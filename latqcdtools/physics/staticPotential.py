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
# Constants for V_Teq0 (from Phys. Rev. D90 (2014) 094503)
#------------------------------------------------------------------------------
# Fit result: [ -91.30436191 1022.25286821  106.70659264]
# Fit error:  [0.53809612 2.51598869 2.58370288]
# chi2/dof:   0.8083127937775374
COULOMB_TERM = -91.30436191  # Coulomb term coefficient [MeV·fm]
STRING_TERM = 1022.25286821  # String tension [MeV/fm]
CONSTANT_TERM = 106.70659264  # Constant term [MeV]

#------------------------------------------------------------------------------
# Static potential functions
#------------------------------------------------------------------------------
def V_Teq0(r: float) -> float:
    """ 
    Zero temperature quark potential in [MeV], takes r in [fm].
    
    The parameters come from a Levenberg-Marquardt fit of the data in 
    Fig 14 of Phys. Rev. D90 (2014) 094503. These values can be obtained
    by running analysistoolbox/hisq_potential/fit_hisq_pot.py.
    
    Args:
        r: Distance in femtometers [fm]
        
    Returns:
        Potential in MeV
    """
    return COULOMB_TERM/r + STRING_TERM*r + CONSTANT_TERM


def fitV_Teq0(r: float, a: float, b: float, c: float) -> float:
    """ 
    Standard Cornell potential fit form.
    
    Args:
        r: Distance in lattice units
        a: Constant term
        b: Coulomb term coefficient
        c: String tension
        
    Returns:
        Potential in lattice units
    """
    return a + b/r + c*r


def fitV_Teq0_oneloop(r: float, a: float, b: float, c: float, d: float) -> float:
    """ 
    Cornell potential with one-loop corrections to Coulomb term.
    
    References:
        - Nucl. Phys. B 129 (1977)
        - Phys. Lett B 92 (1980)
    
    Args:
        r: Distance in lattice units
        a: Constant term
        b: Coulomb term coefficient
        c: String tension
        d: Log term coefficient
        
    Returns:
        Potential in lattice units
    """
    return a + (b + d*np.log(r))/r + c*r


def fitV_Teq0_twoloop(r: float, a: float, b: float, c: float, d: float, e: float) -> float:
    """ 
    Cornell potential with two-loop corrections to Coulomb term.
    
    Reference:
        - Nucl. Phys. B 501 (1997)
    
    Args:
        r: Distance in lattice units
        a: Constant term
        b: Coulomb term coefficient
        c: String tension
        d: Log term coefficient
        e: Log-log term coefficient
        
    Returns:
        Potential in lattice units
    """
    return a + (b + d*np.log(r) + e*np.log(np.log(r)))/r + c*r


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
def compute_sin_terms_kernel(sine_values, improved_coeff, lattice_size, sin_terms_lookup):
    """
    Precompute sin terms for all k1, k2, k3 combinations to avoid redundant calculations.
    
    Args:
        sine_values: Precomputed sine values
        improved_coeff: Coefficient for improved action
        lattice_size: Spatial extension of lattice
        sin_terms_lookup: Output array for precomputed sin terms
    """
    # Get thread position
    idx = cuda.grid(1)
    
    # Calculate 3D indices
    total_indices = lattice_size**3
    if idx >= total_indices:
        return
    
    k3 = idx % lattice_size
    temp = idx // lattice_size
    k2 = temp % lattice_size
    k1 = temp // lattice_size
    
    if not (k1 + k2 + k3 == 0):
        # Calculate sin term once and store it
        sin_term = 0.0
        for k in [k1, k2, k3]:
            sin_k = sine_values[k]
            sin_term += sin_k*sin_k + improved_coeff*sin_k*sin_k*sin_k*sin_k
        
        # Store result in 3D lookup table
        sin_terms_lookup[k1, k2, k3] = 1.0 / (4.0 * sin_term)
    else:
        sin_terms_lookup[k1, k2, k3] = 0.0


@cuda.jit
def calculate_potential_kernel_optimized(points, lattice_size, cosine_values, 
                                        sin_terms_lookup, point_results):
    """
    Optimized CUDA kernel with reduced atomic operations and improved memory access.
    
    Args:
        points: Array of points (x, y, z) to compute, shape (N, 3)
        lattice_size: Spatial extension of lattice
        cosine_values: Precomputed cosine values
        sin_terms_lookup: Precomputed sin terms for all k1,k2,k3 combinations
        point_results: Output array for results (potential, weight) per point
    """
    # Use shared memory for frequently accessed values
    shared_dist_squared = cuda.shared.array(shape=1, dtype=np.int32)
    shared_potential = cuda.shared.array(shape=1, dtype=np.float64)
    
    # Get thread indices
    point_idx = cuda.blockIdx.x
    thread_idx = cuda.threadIdx.x
    
    # Total threads in block
    block_size = cuda.blockDim.x
    
    # Check if this thread's block should process a point
    if point_idx >= points.shape[0]:
        return
    
    # Unpack coordinates (all threads in block access same point)
    x = points[point_idx, 0]
    y = points[point_idx, 1]
    z = points[point_idx, 2]
    
    # Calculate squared distance (only once per block)
    if thread_idx == 0:
        shared_dist_squared[0] = x*x + y*y + z*z
        shared_potential[0] = 0.0
    
    cuda.syncthreads()
    
    # Compute potential with work distributed across threads in block
    # Each thread handles a subset of the k1,k2,k3 combinations
    local_potential = 0.0
    
    # Determine the range each thread processes in the k1,k2,k3 loops
    total_iterations = lattice_size**3
    iterations_per_thread = (total_iterations + block_size - 1) // block_size
    start_idx = thread_idx * iterations_per_thread
    end_idx = min(start_idx + iterations_per_thread, total_iterations)
    
    for idx in range(start_idx, end_idx):
        # Convert linear index to 3D indices
        k3 = idx % lattice_size
        temp = idx // lattice_size
        k2 = temp % lattice_size
        k1 = temp // lattice_size
        
        # Skip origin point
        if not (k1 + k2 + k3 == 0):
            # Index for cosine lookup 
            cos_idx = (k1*x + k2*y + k3*z) % (3*lattice_size*lattice_size)
            cos_term = cosine_values[cos_idx]
            
            # Get precomputed sin term
            sin_term_factor = sin_terms_lookup[k1, k2, k3]
            
            # Add contribution
            local_potential += cos_term * sin_term_factor
    
    # Use shared memory to accumulate results within block
    cuda.atomic.add(shared_potential, 0, local_potential)
    cuda.syncthreads()
    
    # Only first thread in block writes the final result
    if thread_idx == 0:
        # Normalize and store result for this point
        final_potential = shared_potential[0] / lattice_size**3
        dist_squared = shared_dist_squared[0]
        
        # Store results for this point (no atomic operations across points)
        point_results[point_idx, 0] = dist_squared
        point_results[point_idx, 1] = final_potential
        point_results[point_idx, 2] = 1.0  # Weight


@cuda.jit
def calculate_x_potential_kernel_optimized(points, lattice_size, cosine_values,
                                          sin_terms_lookup, point_results):
    """
    Optimized CUDA kernel for x-axis points with reduced atomic operations.
    
    Args:
        points: Array of points (x, 0, 0) to compute, shape (N, 3)
        lattice_size: Spatial extension of lattice
        cosine_values: Precomputed cosine values
        sin_terms_lookup: Precomputed sin terms for all k1,k2,k3 combinations
        point_results: Output array for results (potential, weight) per point
    """
    # Use shared memory for frequently accessed values
    shared_dist_squared = cuda.shared.array(shape=1, dtype=np.int32)
    shared_potential = cuda.shared.array(shape=1, dtype=np.float64)
    
    # Get thread indices
    point_idx = cuda.blockIdx.x
    thread_idx = cuda.threadIdx.x
    
    # Total threads in block
    block_size = cuda.blockDim.x
    
    # Check if this thread's block should process a point
    if point_idx >= points.shape[0]:
        return
    
    # Unpack x coordinate (all threads in block access same point)
    x = points[point_idx, 0]
    
    # Calculate squared distance (only once per block)
    if thread_idx == 0:
        shared_dist_squared[0] = x*x
        shared_potential[0] = 0.0
    
    cuda.syncthreads()
    
    # Compute potential with work distributed across threads in block
    # Each thread handles a subset of the k1,k2,k3 combinations
    local_potential = 0.0
    
    # Determine the range each thread processes in the k1,k2,k3 loops
    total_iterations = lattice_size**3
    iterations_per_thread = (total_iterations + block_size - 1) // block_size
    start_idx = thread_idx * iterations_per_thread
    end_idx = min(start_idx + iterations_per_thread, total_iterations)
    
    for idx in range(start_idx, end_idx):
        # Convert linear index to 3D indices
        k3 = idx % lattice_size
        temp = idx // lattice_size
        k2 = temp % lattice_size
        k1 = temp // lattice_size
        
        # Skip origin point
        if not (k1 + k2 + k3 == 0):
            # Index for cosine lookup (only x component for x-axis points)
            cos_idx = (k1*x) % (3*lattice_size*lattice_size)
            cos_term = cosine_values[cos_idx]
            
            # Get precomputed sin term
            sin_term_factor = sin_terms_lookup[k1, k2, k3]
            
            # Add contribution
            local_potential += cos_term * sin_term_factor
    
    # Use shared memory to accumulate results within block
    cuda.atomic.add(shared_potential, 0, local_potential)
    cuda.syncthreads()
    
    # Only first thread in block writes the final result
    if thread_idx == 0:
        # Normalize and store result for this point
        final_potential = shared_potential[0] / lattice_size**3
        dist_squared = shared_dist_squared[0]
        
        # Store results for this point (no atomic operations across points)
        point_results[point_idx, 0] = dist_squared
        point_results[point_idx, 1] = final_potential
        point_results[point_idx, 2] = 1.0  # Weight


@cuda.jit
def reduce_results_kernel(point_results, potentials, weights, num_points):
    """
    Reduce individual point results into the final potentials and weights arrays.
    
    Args:
        point_results: Array of (dist_squared, potential, weight) for each point
        potentials: Output array for potential values
        weights: Output array for weight counts
        num_points: Number of points
    """
    # Get thread index
    idx = cuda.grid(1)
    
    if idx < num_points:
        dist_squared = int(point_results[idx, 0])
        potential = point_results[idx, 1]
        weight = point_results[idx, 2]
        
        # Update global arrays with atomic operations
        # This is more efficient as we do it once per point instead of in the main kernel
        cuda.atomic.add(potentials, dist_squared, potential)
        cuda.atomic.add(weights, dist_squared, weight)


#------------------------------------------------------------------------------
# CPU fallback implementation
#------------------------------------------------------------------------------
def _cpu_impdist(lattice_size, max_dist_squared, improved_action=True):
    """
    CPU implementation of improved distances calculation.
    Used as fallback when CUDA is not available.
    
    Args:
        lattice_size: Spatial extension of lattice
        max_dist_squared: Maximum squared distance to improve
        improved_action: Whether to use improved action
        
    Returns:
        List of improved distances
    """
    # Set coefficient based on action type
    improved_coeff = 1/3 if improved_action else 0

    @compile
    def compiled_imp_dist():
        """ 
        Ported from code by O. Kaczmarek. 
        """
        improved_distances = []
        k_norm = 2.0 * np.pi / lattice_size
        potentials = [0.0] * 3 * lattice_size**2
        weights = [0] * 3 * lattice_size**2
        
        # Precompute trigonometric values
        cosine_values = [np.cos(i * k_norm) for i in range(3 * lattice_size**2)]
        sine_values = [np.sin(i * k_norm / 2.0) for i in range(lattice_size)]
        
        # First part: compute for all x,y,z points up to Ns/4
        for x in range(int(lattice_size/4) + 1):
            for y in range(int(lattice_size/4) + 1):
                for z in range(int(lattice_size/4) + 1):
                    dist_squared = x**2 + y**2 + z**2
                    if dist_squared > max_dist_squared:
                        continue
                        
                    potential = 0.0
                    for k1 in range(lattice_size):
                        for k2 in range(lattice_size):
                            for k3 in range(lattice_size):
                                if not (k1+k2+k3) == 0:
                                    cos_term = cosine_values[k1*x + k2*y + k3*z]
                                    sin_term = (sine_values[k1]**2 + improved_coeff*sine_values[k1]**4 +
                                               sine_values[k2]**2 + improved_coeff*sine_values[k2]**4 +
                                               sine_values[k3]**2 + improved_coeff*sine_values[k3]**4)
                                    potential += cos_term / (4.0 * sin_term)
                                    
                    potential *= 1.0 / lattice_size**3
                    potentials[dist_squared] += potential
                    weights[dist_squared] += 1
        
        # Second part: compute for x-axis points from Ns/4+1 to Ns/2
        for x in range(int(lattice_size/4) + 1, int(lattice_size/2) + 1):
            dist_squared = x**2
            if dist_squared > max_dist_squared:
                continue
                
            potential = 0.0
            for k1 in range(lattice_size):
                for k2 in range(lattice_size):
                    for k3 in range(lattice_size):
                        if not (k1+k2+k3) == 0:
                            cos_term = cosine_values[k1*x]
                            sin_term = (sine_values[k1]**2 + improved_coeff*sine_values[k1]**4 +
                                       sine_values[k2]**2 + improved_coeff*sine_values[k2]**4 +
                                       sine_values[k3]**2 + improved_coeff*sine_values[k3]**4)
                            potential += cos_term / (4.0 * sin_term)
                            
            potential *= 1.0 / lattice_size**3
            potentials[dist_squared] += potential
            weights[dist_squared] += 1
        
        # Calculate final improved distances
        for i in range(1, int(max_dist_squared) + 1):
            if weights[i] != 0:
                improved_distances.append(
                    1.0 / (4.0 * np.pi * (potentials[i] / weights[i] + 0.22578 / lattice_size))
                )
                
        return improved_distances

    return compiled_imp_dist()

def impdist(lattice_size, max_dist_squared, improved_action=True):
    """
    GPU-accelerated calculation of tree-level improved distances.
    
    Follows equation (3) of 10.1103/PhysRevD.90.074038.
    Falls back to CPU implementation if CUDA is unavailable.

    Args:
        lattice_size: Spatial extension of lattice
        max_dist_squared: Maximum squared distance to improve
        improved_action: Whether to use improved action (default: True)

    Returns:
        List of improved distances

    Raises:
        ValueError: If lattice_size <= 0 or max_dist_squared is too large
    """
    # Input validation
    if not lattice_size > 0:
        logger.TBError(f"Lattice size must be positive, got {lattice_size}")
        
    if max_dist_squared > (lattice_size/2)**2:
        logger.TBError(f"Maximum squared distance {max_dist_squared} exceeds allowed range {(lattice_size/2)**2}")
    
    # Check if CUDA is available
    try:
        if not cuda.is_available():
            logger.warn("CUDA not available, falling back to CPU implementation.")
            return _cpu_impdist(lattice_size, max_dist_squared, improved_action)
    except:
        logger.warn("CUDA not available, falling back to CPU implementation.")
        return _cpu_impdist(lattice_size, max_dist_squared, improved_action)
    
    # Set coefficients based on improved action
    improved_coeff = 1/3 if improved_action else 0
    
    # Prepare data
    k_norm = 2.0 * np.pi / lattice_size
    
    # Pre-calculate cosine and sine values
    cosine_values = np.array([np.cos(i*k_norm) for i in range(3*lattice_size**2)], dtype=np.float64)
    sine_values = np.array([np.sin(i*k_norm/2.0) for i in range(lattice_size)], dtype=np.float64)
    
    # Initialize arrays for potentials and weights
    potentials = np.zeros(3*lattice_size**2, dtype=np.float64)
    weights = np.zeros(3*lattice_size**2, dtype=np.int32)
    
    # Transfer data to device
    cosine_device = cuda.to_device(cosine_values)
    sine_device = cuda.to_device(sine_values)
    potentials_device = cuda.to_device(potentials)
    weights_device = cuda.to_device(weights)
    
    # Get optimal threads per block
    threads_per_block = get_optimal_block_size()
    
    # Create and precompute sin_terms_lookup table to avoid redundant calculations
    sin_terms_lookup = np.zeros((lattice_size, lattice_size, lattice_size), dtype=np.float64)
    sin_terms_device = cuda.to_device(sin_terms_lookup)
    
    # Configure grid for sin terms precomputation
    precompute_threads = threads_per_block
    precompute_blocks = (lattice_size**3 + precompute_threads - 1) // precompute_threads
    
    # Launch kernel to precompute sin terms
    compute_sin_terms_kernel[precompute_blocks, precompute_threads](
        sine_device, improved_coeff, lattice_size, sin_terms_device
    )
    
    # Generate xyz points for first calculation phase
    xyz_points = []
    for x in range(int(lattice_size/4)+1):
        for y in range(int(lattice_size/4)+1):
            for z in range(int(lattice_size/4)+1):
                dist_squared = x**2 + y**2 + z**2
                if dist_squared <= max_dist_squared:
                    xyz_points.append((x, y, z))
    
    if xyz_points:
        # Convert to numpy array and transfer to device
        xyz_points_array = np.array(xyz_points, dtype=np.int32)
        num_xyz_points = len(xyz_points)
        xyz_points_device = cuda.to_device(xyz_points_array)
        
        # Allocate array for individual point results
        xyz_results = np.zeros((num_xyz_points, 3), dtype=np.float64)  # [dist_squared, potential, weight]
        xyz_results_device = cuda.to_device(xyz_results)
        
        # Configure CUDA grid - one block per point, multiple threads per block
        # Each point is processed by an entire block of threads for parallelism
        blocks_per_grid = num_xyz_points
        
        # Launch optimized kernel - one block per point
        calculate_potential_kernel_optimized[blocks_per_grid, threads_per_block](
            xyz_points_device, lattice_size, cosine_device, sin_terms_device, xyz_results_device
        )
        
        # Reduce individual point results to final arrays
        reduce_blocks = (num_xyz_points + threads_per_block - 1) // threads_per_block
        reduce_results_kernel[reduce_blocks, threads_per_block](
            xyz_results_device, potentials_device, weights_device, num_xyz_points
        )
    
    # Generate x-only points for second calculation phase
    x_points = []
    for x in range(int(lattice_size/4)+1, int(lattice_size/2)+1):
        dist_squared = x**2
        if dist_squared <= max_dist_squared:
            x_points.append((x, 0, 0))
    
    if x_points:
        # Convert to numpy array and transfer to device
        x_points_array = np.array(x_points, dtype=np.int32)
        num_x_points = len(x_points)
        x_points_device = cuda.to_device(x_points_array)
        
        # Allocate array for individual point results
        x_results = np.zeros((num_x_points, 3), dtype=np.float64)  # [dist_squared, potential, weight]
        x_results_device = cuda.to_device(x_results)
        
        # Configure CUDA grid - one block per point, multiple threads per block
        blocks_per_grid = num_x_points
        
        # Launch optimized kernel
        calculate_x_potential_kernel_optimized[blocks_per_grid, threads_per_block](
            x_points_device, lattice_size, cosine_device, sin_terms_device, x_results_device
        )
        
        # Reduce individual point results to final arrays
        reduce_blocks = (num_x_points + threads_per_block - 1) // threads_per_block
        reduce_results_kernel[reduce_blocks, threads_per_block](
            x_results_device, potentials_device, weights_device, num_x_points
        )
    
    # Transfer results back to host
    potentials_device.copy_to_host(potentials)
    weights_device.copy_to_host(weights)
    
    # Calculate improved distances
    improved_distances = []
    for i in range(1, int(max_dist_squared)+1):
        if weights[i] != 0:
            improved_distances.append(
                1.0 / (4.0 * np.pi * (potentials[i] / weights[i] + 0.22578/lattice_size))
            )
    
    return improved_distances