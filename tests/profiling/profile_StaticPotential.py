#
# profile_StaticPotential
# 
# K. Ebira 
# 
# Benchmarking script to test the performance scaling of impdist
# with increasing lattice sizes and fixed r2max
#

import matplotlib.pyplot as plt
from latqcdtools.physics.staticPotential import impdist
from latqcdtools.base.plotting import saveFigure
from latqcdtools.base.utilities import timer
import latqcdtools.base.logger as logger
from numba import cuda

def run_benchmarks():
    """
    Run benchmarks for impdist with increasing lattice sizes.
    Compare CPU and GPU performance.
    """
    # Fixed r2max for all tests
    r2max = 20
    
    # Range of lattice sizes to test
    lattice_sizes = [32, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256]
    
    # Check CUDA availability
    cuda_available = False
    try:
        cuda_available = cuda.is_available()
    except:
        pass
    
    logger.info(f"CUDA available: {cuda_available}")
    logger.info(f"Benchmarking impdist with r2max={r2max}")
    
    if cuda_available:
        logger.info(f"{'Lattice Size':<15} {'GPU Time (s)':<15} {'CPU Time (s)':<15} {'Speedup':<10}")
    else:
        logger.info(f"{'Lattice Size':<15} {'Time (s)':<15}")
    
    results = []
    clock   = timer()
    
    for ns in lattice_sizes:
        # Skip if invalid parameters
        if r2max > (ns/2)**2:
            logger.info(f"Skipping Ns={ns} because r2max={r2max} exceeds maximum allowed ({(ns/2)**2})")
            continue
        
        if cuda_available:
            # For CPU timing, temporarily patch cuda.is_available
            original_is_available = cuda.is_available
            cuda.is_available = lambda: False
            
            # Time CPU implementation
            clock.resetTimer()
            _ = impdist(ns, r2max)
            cpu_time = clock.printTiming('CPU')
            
            # Restore CUDA availability check
            cuda.is_available = original_is_available
            
            # Time GPU implementation
            clock.resetTimer()
            _ = impdist(ns, r2max)
            gpu_time = clock.printTiming('GPU')
            
            # Calculate speedup
            speedup = cpu_time / gpu_time
            
            logger.info(f"{ns:<15} {gpu_time:<15.6f} {cpu_time:<15.6f} {speedup:<10.2f}x")
            results.append((ns, gpu_time, cpu_time, speedup))
        else:
            # Just time the CPU implementation (only option)
            clock.resetTimer()
            _ = impdist(ns, r2max)
            elapsed = clock.printTiming('CPU')
            
            logger.info(f"{ns:<15} {elapsed:<15.6f}")
            results.append((ns, None, elapsed, None))
    
    return results

def plot_results(results):
    """
    Plot benchmark results.
    
    Args:
        results: List of tuples (lattice_size, gpu_time, cpu_time, speedup)
    """
    lattice_sizes = [r[0] for r in results]
    gpu_times = [r[1] for r in results]
    cpu_times = [r[2] for r in results]
    speedups = [r[3] for r in results]
    
    has_gpu_data = gpu_times[0] is not None
    
    plt.figure(figsize=(12, 10))
    
    # Plot 1: Runtime vs Lattice Size
    plt.subplot(2, 1, 1)
    plt.plot(lattice_sizes, cpu_times, 'o-', label='CPU', linewidth=2, markersize=8)
    
    if has_gpu_data:
        plt.plot(lattice_sizes, gpu_times, 's-', label='GPU', linewidth=2, markersize=8)
    
    plt.xlabel('Lattice Size (Ns)')
    plt.ylabel('Runtime (seconds)')
    plt.title('Strong Scaling: Runtime vs Lattice Size (r2max=20)')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    
    # Plot 2: Speedup vs Lattice Size (if GPU data available)
    if has_gpu_data:
        plt.subplot(2, 1, 2)
        plt.plot(lattice_sizes, speedups, 'D-', color='green', linewidth=2, markersize=8)
        plt.axhline(y=1.0, color='r', linestyle='--', alpha=0.7)
        plt.xlabel('Lattice Size (Ns)')
        plt.ylabel('Speedup (CPU time / GPU time)')
        plt.title('GPU Speedup vs Lattice Size')
        plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    saveFigure('static_potential_scaling.png', dpi=300)
    plt.show()

if __name__ == '__main__':
    results = run_benchmarks()
    plot_results(results)
