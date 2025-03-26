# staticPotential.py
# 
# D. Clarke, K. Ebira
# 
# Some methods related to the static potential in lattice QCD.
# Optimized GPU implementation for improved scaling with lattice size.

import numpy as np
from numba import cuda
import latqcdtools.base.logger as logger
from latqcdtools.base.speedify import numbaON, compile

numbaON()

#------------------------------------------------------------------------------
# Static potential functions
#------------------------------------------------------------------------------
def V_Teq0(r) -> float:
    """ 
    Zero temperature quark potential in [MeV], takes r in [fm]. The parameters a, b, and c come from
    a Levenberg-Marquardt fit of the data in Fig 14 of Phys. Rev. D90 (2014) 094503. These numbers can
    be obtained again by running analysistoolbox/hisq_potential/fit_hisq_pot.py. 
    """
    #    result:  [ -91.30436191 1022.25286821  106.70659264]
    #     error:  [0.53809612 2.51598869 2.58370288]
    #  chi2/dof:  0.8083127937775374
    a =  -91.30436191
    b = 1022.25286821
    c =  106.70659264
    return a/r + b*r + c



def fitV_Teq0(r, a, b, c) -> float:
    """ 
    Fit form of standard Cornell potential. Fit to be done in lattice units.
    """
    return a + b/r + c*r


def fitV_Teq0_oneloop(r,a,b,c,d) -> float:
    """ 
    Including one-loop corrections to Coulomb. See Nucl. Phys. B 129 (1977) and Phys. Lett B 92 (1980).
    Fit to be done in lattice units.
    """
    return a + ( b + d*np.log(r) )/r + c*r


def fitV_Teq0_twoloop(r,a,b,c,d,e) -> float:
    """ 
    Including two-loop corrections to Coulomb. See Nucl. Phys. B 501 (1997). Fit to be done in lattice units.
    """
    return a + ( b + d*np.log(r) + e*np.log(np.log(r)) )/r + c*r


#------------------------------------------------------------------------------
# CUDA utility functions
#------------------------------------------------------------------------------
def get_optimal_block_size():
    """
    Returns an optimal block size based on the current CUDA device.
    Defaults to 256 if device information cannot be obtained.
    
    Returns:
        int: Optimal threads per block
    """
    try:
        device = cuda.get_current_device()
        # Max threads per block is device.MAX_THREADS_PER_BLOCK, but usually better to use less
        return min(256, device.MAX_THREADS_PER_BLOCK)
    except:
        return 256


#------------------------------------------------------------------------------
# CUDA kernels
#------------------------------------------------------------------------------
@cuda.jit
def compute_sine_terms_kernel(sinf, cw, Ns, sine_terms_lookup):
    """
    Precompute sin terms for all k1, k2, k3 combinations to improve performance.
    
    Args:
        sinf: Precomputed sine values
        cw: Coefficient for improved action
        Ns: Spatial extension of lattice
        sine_terms_lookup: Output array for precomputed sine terms
    """
    # Get thread position
    idx = cuda.grid(1)
    
    # Calculate 3D indices
    total_indices = Ns**3
    if idx >= total_indices:
        return
    
    k3 = idx % Ns
    temp = idx // Ns
    k2 = temp % Ns
    k1 = temp // Ns
    
    if not (k1 + k2 + k3 == 0):
        # Calculate sine term once and store it
        r2 = 0.0
        for k in [k1, k2, k3]:
            sin_k = sinf[k]
            r2 += sin_k*sin_k + cw*sin_k*sin_k*sin_k*sin_k
        
        # Store result in 3D lookup table
        sine_terms_lookup[k1, k2, k3] = 1.0 / (4.0 * r2)
    else:
        sine_terms_lookup[k1, k2, k3] = 0.0

@cuda.jit
def calculate_potential_kernel(points, Ns, cosf, sine_terms_lookup, point_results):
    """
    CUDA kernel for calculating lattice QCD potentials for a set of spatial points.
    
    This kernel calculates the static potential for each input point. Each block 
    processes one spatial point, with the workload distributed across threads within the block.
    
    Args:
        points: Array of spatial coordinates (x,y,z) with shape (N,3)
        Ns: Size of the lattice in each dimension
        cosf: Precomputed cosine values for all possible indices
        sine_terms_lookup: Precomputed sine terms for all (k1,k2,k3) combinations
        point_results: Output array to store results with shape (N,3)
                       Each row contains [sq, pot, weight]
    """
    # Allocate shared memory for values accessed by all threads in block
    shared_sq = cuda.shared.array(shape=1, dtype=np.int32)
    shared_pot = cuda.shared.array(shape=1, dtype=np.float64)
    
    # Get the current thread's block and thread indices
    point_idx = cuda.blockIdx.x
    thread_idx = cuda.threadIdx.x
    
    # Get number of threads in this block
    block_size = cuda.blockDim.x
    
    # Exit if this block doesn't correspond to a valid point
    if point_idx >= points.shape[0]:
        return
    
    # Get coordinates of the point for this block
    x = points[point_idx, 0]
    y = points[point_idx, 1]
    z = points[point_idx, 2]
    
    # Thread 0 initializes shared memory values
    if thread_idx == 0:
        shared_sq[0] = x*x + y*y + z*z
        shared_pot[0] = 0.0
    
    # Wait for thread 0 to complete initialization
    cuda.syncthreads()
    
    # Each thread will calculate a partial sum of the potential
    local_pot = 0.0
    
    # Divide the workload among threads in the block
    total_iterations = Ns**3
    iterations_per_thread = (total_iterations + block_size - 1) // block_size
    start_idx = thread_idx * iterations_per_thread
    end_idx = min(start_idx + iterations_per_thread, total_iterations)
    
    # Calculate this thread's portion of the potential sum
    for idx in range(start_idx, end_idx):
        # Convert linear index to 3D lattice indices
        k3 = idx % Ns
        temp = idx // Ns
        k2 = temp % Ns
        k1 = temp // Ns
        
        # Skip the origin point in momentum space
        if not (k1 + k2 + k3 == 0):
            # Calculate index for cosine lookup
            cos_idx = (k1*x + k2*y + k3*z) % (3*Ns*Ns)
            r1 = cosf[cos_idx]
            
            # Get the precomputed sine term for these k-indices
            r2_factor = sine_terms_lookup[k1, k2, k3]
            
            # Add this term's contribution to the local sum
            local_pot += r1 * r2_factor
    
    # Add local result to the block's shared sum
    cuda.atomic.add(shared_pot, 0, local_pot)
    
    # Wait for all threads to complete their additions
    cuda.syncthreads()
    
    # Thread 0 finalizes and stores the results
    if thread_idx == 0:
        # Normalize the potential by the lattice volume
        final_pot = shared_pot[0] / Ns**3
        sq = shared_sq[0]
        
        # Store results for this point
        point_results[point_idx, 0] = sq
        point_results[point_idx, 1] = final_pot
        point_results[point_idx, 2] = 1.0  # Weight for later averaging


@cuda.jit
def calculate_x_potential_kernel(points, Ns, cosf,
                                 sine_terms_lookup, point_results):
    """
    CUDA kernel for calculating potentials for points on the x-axis.
    
    This kernel computes the static potential specifically for points along the 
    x-axis (where y=0, z=0). 
    
    Each block processes one spatial point with work distributed across threads.
    
    Args:
        points: Array of x-axis coordinates with shape (N,3), where the y and z
                components are expected to be zero
        Ns: Size of the lattice in each dimension
        cosf: Precomputed cosine values for all possible indices
        sine_terms_lookup: Precomputed sine terms for all (k1,k2,k3) combinations
        point_results: Output array to store results with shape (N,3)
                       Each row contains [sq, pot, weight]
    """
    # Allocate shared memory for values accessed by all threads in block
    shared_sq = cuda.shared.array(shape=1, dtype=np.int32)
    shared_pot = cuda.shared.array(shape=1, dtype=np.float64)
    
    # Get the current thread's block and thread indices
    point_idx = cuda.blockIdx.x
    thread_idx = cuda.threadIdx.x
    
    # Get number of threads in this block
    block_size = cuda.blockDim.x
    
    # Exit if this block doesn't correspond to a valid point
    if point_idx >= points.shape[0]:
        return
    
    # Get x-coordinate of the point for this block (y and z are zero)
    x = points[point_idx, 0]
    
    # Thread 0 initializes shared memory values
    if thread_idx == 0:
        shared_sq[0] = x*x  # Simplified: sq = x**2 for x-axis points
        shared_pot[0] = 0.0
    
    # Wait for thread 0 to complete initialization
    cuda.syncthreads()
    
    # Each thread will calculate a partial sum of the potential
    local_pot = 0.0
    
    # Divide the workload among threads in the block
    total_iterations = Ns**3
    iterations_per_thread = (total_iterations + block_size - 1) // block_size
    start_idx = thread_idx * iterations_per_thread
    end_idx = min(start_idx + iterations_per_thread, total_iterations)
    
    # Calculate this thread's portion of the potential sum
    for idx in range(start_idx, end_idx):
        # Convert linear index to 3D lattice indices
        k3 = idx % Ns
        temp = idx // Ns
        k2 = temp % Ns
        k1 = temp // Ns
        
        # Skip the origin point in momentum space
        if not (k1 + k2 + k3 == 0):
            # Simplified cosine lookup for x-axis points (y=0, z=0)
            cos_idx = (k1*x) % (3*Ns*Ns)
            r1 = cosf[cos_idx]
            
            # Get the precomputed sine term for these k-indices
            r2_factor = sine_terms_lookup[k1, k2, k3]
            
            # Add this term's contribution to the local sum
            local_pot += r1 * r2_factor
    
    # Add local result to the block's shared sum
    cuda.atomic.add(shared_pot, 0, local_pot)
    
    # Wait for all threads to complete their additions
    cuda.syncthreads()
    
    # Thread 0 finalizes and stores the results
    if thread_idx == 0:
        # Normalize the potential by the lattice volume
        final_pot = shared_pot[0] / Ns**3
        sq = shared_sq[0]
        
        # Store results for this point
        point_results[point_idx, 0] = sq
        point_results[point_idx, 1] = final_pot
        point_results[point_idx, 2] = 1.0  # Weight for later averaging


@cuda.jit
def reduce_results_kernel(point_results, pots, weight, num_points):
    """
    CUDA kernel for aggregating individual point calculations into final result arrays.
    
    This kernel takes the individual point results from
    previous calculations and combines them by their squared distance, summing both
    the potential values and weights for subsequent averaging.
    
    The reduction process organizes results by distance rather than by spatial 
    coordinates, allowing the static potential to be properly represented as a
    function of separation distance.
    
    Args:
        point_results: Array of results from individual points with shape (N,3),
                      where each row contains [sq, pot, weight]
        pots: Output array indexed by squared distance, accumulating potential values
        weight: Output array indexed by squared distance, accumulating weights
        num_points: Total number of points to process
    """
    # Get the unique global thread index
    idx = cuda.grid(1)
    
    # Only process if this thread corresponds to a valid point
    if idx < num_points:
        # Extract data for this point
        sq = int(point_results[idx, 0])  # Convert to integer for array indexing
        pot = point_results[idx, 1]      # Potential value for this point
        w = point_results[idx, 2]        # Weight for this point (typically 1.0)
        
        # Atomically add this point's values to the appropriate distance bin
        # Atomic operations are necessary to handle the case where multiple threads
        # attempt to update the same distance bin simultaneously
        cuda.atomic.add(pots, sq, pot)
        cuda.atomic.add(weight, sq, w)


#------------------------------------------------------------------------------
# CPU fallback implementation
#------------------------------------------------------------------------------
def _cpu_impdist(Ns, r2max, improvedAction=True):
    """
    CPU implementation of improved distances calculation.
    Used as fallback when CUDA is not available.
    
    Args:
        Ns: Spatial extension of lattice
        r2max: Maximum squared distance to improve
        improvedAction: Whether to use improved action
        
    Returns:
        List of improved distances
    """
    # Set coefficient based on action type
    cw = 1/3 if improvedAction else 0

    @compile
    def compiled_impdist():
        """ 
        Ported from code by O. Kaczmarek. 
        """
        rimp = []
        kn = 2.0 * np.pi / Ns
        pots = [0.0] * 3 * Ns**2
        weight = [0] * 3 * Ns**2
        
        # Precompute trigonometric values
        cosf = [np.cos(i * kn) for i in range(3 * Ns**2)]
        sinf = [np.sin(i * kn / 2.0) for i in range(Ns)]
        
        # First compute for all x,y,z points up to Ns/4
        for x in range(int(Ns/4) + 1):
            for y in range(int(Ns/4) + 1):
                for z in range(int(Ns/4) + 1):
                    sq = x**2 + y**2 + z**2
                    if sq > r2max:
                        continue
                        
                    pot = 0.0
                    for k1 in range(Ns):
                        for k2 in range(Ns):
                            for k3 in range(Ns):
                                if not (k1+k2+k3) == 0:
                                    r1 = cosf[k1*x + k2*y + k3*z]
                                    r2 = (sinf[k1]**2 + cw*sinf[k1]**4 +
                                          sinf[k2]**2 + cw*sinf[k2]**4 +
                                          sinf[k3]**2 + cw*sinf[k3]**4)
                                    pot += r1 / (4.0 * r2)
                                    
                    pot *= 1.0 / Ns**3
                    pots[sq] += pot
                    weight[sq] += 1
        
        # Then compute for x-axis points from Ns/4+1 to Ns/2
        for x in range(int(Ns/4) + 1, int(Ns/2) + 1):
            sq = x**2
            if sq > r2max:
                continue
                
            pot = 0.0
            for k1 in range(Ns):
                for k2 in range(Ns):
                    for k3 in range(Ns):
                        if not (k1+k2+k3) == 0:
                            r1 = cosf[k1*x]
                            r2 = (sinf[k1]**2 + cw*sinf[k1]**4 +
                                  sinf[k2]**2 + cw*sinf[k2]**4 +
                                  sinf[k3]**2 + cw*sinf[k3]**4)
                            pot += r1 / (4.0 * r2)
                            
            pot *= 1.0 / Ns**3
            pots[sq] += pot
            weight[sq] += 1
        
        # Calculate final improved distances
        for i in range(1, int(r2max) + 1):
            if weight[i] != 0:
                rimp.append(
                    1.0 / (4.0 * np.pi * (pots[i] / weight[i] + 0.22578 / Ns))
                )
                
        return rimp

    return compiled_impdist()

def impdist(Ns, r2max, improvedAction=True):
    """
    GPU-accelerated calculation of tree-level improved distances.
    
    Follows equation (3) of 10.1103/PhysRevD.90.074038.
    Falls back to CPU implementation if CUDA is unavailable.

    Args:
        Ns: Spatial extension of lattice
        r2max: Maximum squared distance to improve
        improvedAction: Whether to use improved action (default: True)

    Returns:
        rimp: List of improved distances

    Raises:
        ValueError: If Ns <= 0 or r2max is too large
    """
    # Input validation
    if not Ns > 0:
        logger.TBError("Need Ns>0")
        
    if r2max > (Ns/2)**2:
        logger.TBError("r2max is too large.")
    
    # Check if CUDA is available
    try:
        if not cuda.is_available():
            logger.warn("CUDA not available, falling back to CPU implementation.")
            return _cpu_impdist(Ns, r2max, improvedAction)
    except:
        logger.warn("CUDA not available, falling back to CPU implementation.")
        return _cpu_impdist(Ns, r2max, improvedAction)

    # Set coefficients based on improved action
    cw = 1/3 if improvedAction else 0
    
    # Prepare data
    kn = 2.0 * np.pi / Ns
    
    # Pre-calculate cosine and sine values
    cosf = np.array([np.cos(i*kn) for i in range(3*Ns**2)], dtype=np.float64)
    sinf = np.array([np.sin(i*kn/2.0) for i in range(Ns)], dtype=np.float64)
    
    # Initialize arrays for potentials and weights
    pots = np.zeros(3*Ns**2, dtype=np.float64)
    weight = np.zeros(3*Ns**2, dtype=np.int32)
    
    # Transfer data to device
    cosf_device = cuda.to_device(cosf)
    sinf_device = cuda.to_device(sinf)
    pots_device = cuda.to_device(pots)
    weight_device = cuda.to_device(weight)
    
    # Get optimal threads per block
    threads_per_block = get_optimal_block_size()
    
    # Create and precompute sine_terms_lookup table to avoid redundant calculations
    sine_terms_lookup = np.zeros((Ns, Ns, Ns), dtype=np.float64)
    sine_terms_device = cuda.to_device(sine_terms_lookup)
    
    # Configure grid for sine terms precomputation
    precompute_threads = threads_per_block
    precompute_blocks = (Ns**3 + precompute_threads - 1) // precompute_threads
    
    # Launch kernel to precompute sine terms
    try:
        compute_sine_terms_kernel[precompute_blocks, precompute_threads](
            sinf_device, cw, Ns, sine_terms_device
        )
    except Exception as e:
        logger.warn('Kernel launch failed. Exception:')
        logger.warn(e)
        logger.warn("CUDA not available, falling back to CPU implementation.")
        return _cpu_impdist(Ns, r2max, improvedAction)

    # Generate xyz points for first calculation phase
    xyz_points = []
    for x in range(int(Ns/4)+1):
        for y in range(int(Ns/4)+1):
            for z in range(int(Ns/4)+1):
                sq = x**2 + y**2 + z**2
                if sq <= r2max:
                    xyz_points.append((x, y, z))
    
    if xyz_points:
        # Convert to numpy array and transfer to device
        xyz_points_array = np.array(xyz_points, dtype=np.int32)
        num_xyz_points = len(xyz_points)
        xyz_points_device = cuda.to_device(xyz_points_array)
        
        # Allocate array for individual point results
        xyz_results = np.zeros((num_xyz_points, 3), dtype=np.float64)  # [sq, pot, weight]
        xyz_results_device = cuda.to_device(xyz_results)
        
        # Configure CUDA grid - one block per point, multiple threads per block
        # Each point is processed by an entire block of threads for parallelism
        blocks_per_grid = num_xyz_points
        
        # Launch optimized kernel - one block per point
        calculate_potential_kernel[blocks_per_grid, threads_per_block](
            xyz_points_device, Ns, cosf_device, sine_terms_device, xyz_results_device
        )
        
        # Reduce individual point results to final arrays
        reduce_blocks = (num_xyz_points + threads_per_block - 1) // threads_per_block
        reduce_results_kernel[reduce_blocks, threads_per_block](
            xyz_results_device, pots_device, weight_device, num_xyz_points
        )
    
    # Generate x-only points for second calculation phase
    x_points = []
    for x in range(int(Ns/4)+1, int(Ns/2)+1):
        sq = x**2
        if sq <= r2max:
            x_points.append((x, 0, 0))
    
    if x_points:
        # Convert to numpy array and transfer to device
        x_points_array = np.array(x_points, dtype=np.int32)
        num_x_points = len(x_points)
        x_points_device = cuda.to_device(x_points_array)
        
        # Allocate array for individual point results
        x_results = np.zeros((num_x_points, 3), dtype=np.float64)  # [sq, pot, weight]
        x_results_device = cuda.to_device(x_results)
        
        # Configure CUDA grid - one block per point, multiple threads per block
        blocks_per_grid = num_x_points
        
        # Launch optimized kernel
        calculate_x_potential_kernel[blocks_per_grid, threads_per_block](
            x_points_device, Ns, cosf_device, sine_terms_device, x_results_device
        )
        
        # Reduce individual point results to final arrays
        reduce_blocks = (num_x_points + threads_per_block - 1) // threads_per_block
        reduce_results_kernel[reduce_blocks, threads_per_block](
            x_results_device, pots_device, weight_device, num_x_points
        )
    
    # Transfer results back to host
    pots_device.copy_to_host(pots)
    weight_device.copy_to_host(weight)
    
    # Calculate improved distances
    rimp = []
    for i in range(1, int(r2max)+1):
        if weight[i] != 0:
            rimp.append(
                1.0 / (4.0 * np.pi * (pots[i] / weight[i] + 0.22578/Ns))
            )
    
    return rimp