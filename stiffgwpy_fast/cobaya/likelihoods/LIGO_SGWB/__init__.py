"""Optional LVK stochastic-background cross-correlation likelihood.

中文：按需导出 LVK O1/O2/O3 SGWB 交叉相关似然类。
"""

__all__ = ['LVK_SGWB_CC']


def __getattr__(name):
    if name != 'LVK_SGWB_CC':
        raise AttributeError("module %r has no attribute %r" % (__name__, name))
    from .LVK_SGWB_CC import LVK_SGWB_CC

    globals()[name] = LVK_SGWB_CC
    return LVK_SGWB_CC
