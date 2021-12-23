import numpy as np
from SemiBin.cluster import cal_kl

def slow_cal_kl(m, v):
    """
    Calculate KL divergence
    """
    m = np.clip(m, 1e-6, None)
    v = np.clip(v, 1.0, None)
    m1 = m.reshape(1, len(m))
    m2 = m.reshape(len(m), 1)

    v1 = v.reshape(1, len(v))
    v2 = v.reshape(len(v), 1)

    value = np.log(np.sqrt(v1 / v2)) +  np.divide(np.square(m1 - m2) + v2, 2 * v1) - 0.5

    return np.clip(value, 1e-6, 1 - 1e-6)

def test_cal_kl():
    for _ in range(32):
        m = np.random.rand(128)
        v = np.random.rand(128)
        assert np.allclose(slow_cal_kl(m, v), cal_kl(m ,v, use_ne='yes'))
        assert np.allclose(slow_cal_kl(m, v), cal_kl(m ,v, use_ne='no'))
        assert np.allclose(slow_cal_kl(m, v), cal_kl(m ,v, use_ne='auto'))
