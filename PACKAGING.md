# stiffGWpy 发布说明

本项目按照 Python Packaging User Guide 的标准流程构建。`docs/`、测试文件、验证脚本、
CI 配置和历史研究型配置不属于 PyPI 运行时发布内容；运行所需的 Python 模块、Cobaya
适配器和必要数据文件会保留在 wheel 中。

## 用户安装

发布后，用户可以直接安装：

```bash
python -m pip install stiffgwpy
```

默认调用 `LCDM_SG.SGWB_iter()` 使用 `fast` 引擎的 plain-grid 档位。需要更高精度时，
显式使用 `accuracy_mode="production"`；需要独立精度锚点时使用
`engine="reference"`；需要原始回归路径时使用 `engine="lsoda"`。

Cobaya 集成是可选依赖：

```bash
python -m pip install "stiffgwpy[cobaya]"
```

## 发布前检查

在包含 `pyproject.toml` 的项目根目录执行：

```bash
python -m pip install --upgrade build twine
Remove-Item -Recurse -Force build, dist -ErrorAction SilentlyContinue
python -m build
python scripts/verify_distribution.py dist
```

校验脚本会确认同时生成 wheel 和 sdist，并确认 `docs/`、`tests/`、`scripts/`、`.github/`
以及 `mcmc_compare.yaml` 没有进入发布归档。

## TestPyPI 验证

首次发布或较大改动时，先上传到 TestPyPI：

```bash
python -m twine upload --repository testpypi dist/*
python -m pip install --index-url https://test.pypi.org/simple/ --no-deps stiffgwpy
```

上传需要 TestPyPI API token。不要把 token 写入仓库或命令历史。

仓库根目录的 `.pypirc` 属于本地凭据文件，已加入 `.gitignore`，只可用于本地发布，
不能提交到 GitHub；发布完成后应继续保存在本机安全位置。

## 正式 PyPI 发布

确认 TestPyPI 安装和导入成功后，再执行：

```bash
python -m twine upload dist/*
python -m pip install --upgrade stiffgwpy
```

每次正式发布前必须递增 `pyproject.toml` 的 `project.version`，并同步更新
`CHANGELOG.md`、`README.md` 和 `README_zh.md`。当前项目发布邮箱为
`2966684515@qq.com`；PyPI 登录使用 API token，不使用邮箱密码。
