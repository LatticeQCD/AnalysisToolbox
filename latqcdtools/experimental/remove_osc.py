import numpy as np

# See https://arxiv.org/abs/1411.3018

def remove_osc_av(corr, xdata):
    G_no = []
    G_o = []
    m_no = []
    m_o = []
    A_no = []
    A_o = []
    xdata_split = []
    for i in range(len(xdata) - 4):
        A = corr[i + 1]**2 - corr[i + 2] * corr[i]
        B = corr[i + 3] * corr[i] - corr[i + 2] * corr[i + 1]
        C = corr[i + 2]**2 - corr[i + 3] * corr[i + 1]
        x_p = B / (2 * A) + np.sqrt(B**2 - 4 * A * C) / (2 * np.abs(A))
        x_m = -B / (2 * A) + np.sqrt(B**2 - 4 * A * C) / (2 * np.abs(A))
        if x_m < 0:
            x_m = np.nan
        if x_p < 0:
            x_p = np.nan
        m_p = -np.log(x_p)
        m_m = -np.log(x_m)
        m_o.append(m_p)
        m_no.append(m_m)
        G_no.append((corr[i + 1] + corr[i] * x_p) / (x_m + x_p))
        G_o.append(abs((-1)**xdata[i] * (corr[i + 1] - corr[i] * x_m) / (x_m + x_p)))
        A_no.append(G_no[-1] / np.exp(-m_no[-1]*xdata[i]))
        A_o.append(G_o[-1] / np.exp(-m_o[-1]*xdata[i]))
        xdata_split.append(xdata[i])

    return xdata_split, G_no, G_o, m_no, m_o, A_no, A_o

def remove_osc(data, xdata):
    corr = np.array([np.mean(i) for i in data])
    return remove_osc_av(corr, xdata)