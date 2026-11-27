# 仓库脚本

- `verify_distribution.py` 检查 wheel 和源码包不包含仓库专用文件。
- `smoke_installed_wheel.py` 在隔离目录安装 wheel，并检查导入和打包资源。
- `_smoke_probe.py` 运行 CI 使用的 Cobaya 适配器冒烟检查。
- `_resource_budget.py` 及三个核心求解器测试直接导入的数值参考辅助模块作为测试支撑保留，不会进入安装包。

批量扫描、候选原型、绘图、后验工作流和独立基准驱动均保存在本地研究归档。
