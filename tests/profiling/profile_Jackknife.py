# 
# profileJackknife.py                                                               
# 
# D. Clarke, H. Dick 
# 
# See how fast or slow different implementations of the jackknife are 
# 
from latqcdtools.legacy import jackknife as jackknifeLeg
from latqcdtools.statistics.jackknife import jackknife
from latqcdtools.base.speedify import DEFAULTTHREADS
import numpy as np
from latqcdtools.base.utilities import timer
import multiprocessing as mp
import latqcdtools.base.logger as logger

# jackknife estimator for a simple estimator f
# f - function that takes  data array and returns its estimated value (can be an ndarray)
# data - The data array with the shape [samples, ...]
# blocks - the number of jackknife blocks
def jackknife2(f, data, blocks=20):
    data = np.asarray(data)
    block_id = np.linspace(0, blocks, len(data),
                    endpoint=False).astype(np.int32)
    J = [f(data[block_id != i]) for i in range(blocks)]
    # slower version:
    #J = [f(data[np.concatenate([np.arange(len(data)*i//blocks),
    #                 np.arange(len(data)*(i+1)//blocks, len(data))])]) for i in range(blocks)]
    return blocks*f(data) - (blocks-1)*np.mean(J, axis=0),\
           np.std(J, axis=0)*(blocks-1)**.5

# jackknife estimator for a simple estimator f
# f - function that takes  data array and returns its estimated value (can be an ndarray)
# data - The data array with the shape [samples, ...]
# blocks - the number of jackknife blocks
def jackknife3(f, data, blocks=20):
    data = np.asarray(data)
    block_id = np.linspace(0, blocks, len(data),
                    endpoint=False).astype(np.int32)
    # sending the data to the processes (copying it) is very tedious and takes about as much time as the function call itself.
    # this is only worth it for very slow functions f (i.e. if f is O(n^2))
    J = mp.Pool(mp.cpu_count()).map(f, (data[block_id != i] for i in range(blocks)))
    return blocks*f(data) - (blocks-1)*np.mean(J, axis=0),\
           np.std(J, axis=0)*(blocks-1)**.5

data = np.random.random(10000000)

timey = timer()
logger.info('simple')
logger.info(jackknife2(np.std, data, 20))
logger.info(jackknife2(np.median, data, 20))
logger.info(jackknife2(lambda x: 1 / np.mean(x), data, 20))
timey.printTiming()

logger.info('multiprocess')
logger.info(jackknife3(np.std, data, 20))
logger.info(jackknife3(np.median, data, 20))
def inv_mean(x):
    return 1 / np.mean(x)
logger.info(jackknife3(inv_mean, data, 20))
timey.printTiming()

logger.info('legacy')
logger.info(jackknifeLeg(np.std, data, 20, conf_axis=0))
logger.info(jackknifeLeg(np.median, data, 20, conf_axis=0))
def inv_mean(x):
    return 1 / np.mean(x)
logger.info(jackknifeLeg(inv_mean, data, 20, conf_axis=0))
timey.printTiming()

logger.info('current, 1-proc')
logger.info(jackknife(np.std, data, 20, conf_axis=0 ))
logger.info(jackknife(np.median, data, 20, conf_axis=0))
def inv_mean(x):
    return 1 / np.mean(x)
logger.info(jackknife(inv_mean, data, 20, conf_axis=0))
timey.printTiming()

logger.info('current parallel')
logger.info(jackknife(np.std, data, 20, conf_axis=0, nproc=DEFAULTTHREADS))
logger.info(jackknife(np.median, data, 20, conf_axis=0, nproc=DEFAULTTHREADS))
def inv_mean(x):
    return 1 / np.mean(x)
logger.info(jackknife(inv_mean, data, 20, conf_axis=0, nproc=DEFAULTTHREADS))
timey.printTiming()
