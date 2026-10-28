"""Optional IPTA and NANOGrav pulsar-timing likelihood adapters.

中文：按需导出 IPTA 与 NANOGrav 似然类；纯数值 helper 不要求安装 Cobaya。
数据路径和参数默认值由同目录 YAML 提供。
"""

__all__ = ['IPTA', 'NANOGrav']


def __getattr__(name):
    if name == 'IPTA':
        from .IPTA import IPTA

        globals()[name] = IPTA
        return IPTA
    if name == 'NANOGrav':
        from .NANOGrav import NANOGrav

        globals()[name] = NANOGrav
        return NANOGrav
    raise AttributeError("module %r has no attribute %r" % (__name__, name))
