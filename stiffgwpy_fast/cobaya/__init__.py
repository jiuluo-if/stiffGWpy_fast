"""Optional Cobaya theory and likelihood adapters for ``stiffgwpy_fast``.

中文：按需加载 Cobaya 理论适配器；纯数值 likelihood helper 不依赖可选 Cobaya 包。
配置契约见同目录 README 与 `docs/cobaya.md`。
"""

from importlib import import_module

__all__ = ["stiffGW"]


def __getattr__(name):
    if name != "stiffGW":
        raise AttributeError("module %r has no attribute %r" % (__name__, name))
    module = import_module(".stiffGW", __name__)
    adapter = module.stiffGW
    globals()[name] = adapter
    return adapter
