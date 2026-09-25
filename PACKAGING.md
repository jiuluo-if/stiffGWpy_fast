# stiffGWpy_fast 发布说明

本项目使用 Python Packaging User Guide 所述流程构建。发布包包含运行所需的 Python
模块、Cobaya 适配器和数据文件；`docs/`、测试、验证脚本、CI 配置及研究专用配置不属于
运行时包内容。

## 用户安装

发布后，用户可以直接安装：

```bash
python -m pip install stiffgwpy_fast
```

调用 `LCDM_SG.SGWB_iter()` 时默认使用唯一正式的 `fast` 用户档位
（goal-kink-hybrid）。旧的 `plain-grid`、`production` 和 `ultra-fast` 名称不代表其他
正式档位；高层 API 会将这些兼容名称映射到 `fast`。需要独立精度参照时使用
`engine="reference"`；需要原始回归路径时使用 `engine="lsoda"`。`production` 等内部
验证配置不应描述为面向用户的精度档位。

Cobaya 集成是可选依赖：

```bash
python -m pip install "stiffgwpy_fast[cobaya]"
```

## 发布前检查

在包含 `pyproject.toml` 的项目根目录执行：

```bash
python -m pip install --upgrade build twine
python -c "import shutil; shutil.rmtree('build', ignore_errors=True); shutil.rmtree('dist', ignore_errors=True)"
python -m build
python scripts/verify_distribution.py dist
```

清理旧的 `build/` 与 `dist/` 后再构建，避免旧归档干扰检查。校验脚本要求目录中恰有一个
wheel 和一个 sdist，并检查 `docs/`、`tests/`、`scripts/`、`.github/` 及
`mcmc_compare.yaml` 未进入发布归档，同时确认 wheel 含有必要的包文件。

## TestPyPI 验证

首次发布或较大改动时，先上传到 TestPyPI：

```bash
python -m twine upload --repository testpypi dist/*
python -m pip install --index-url https://test.pypi.org/simple/ --no-deps stiffgwpy_fast
```

上传需要 TestPyPI API token。不要把 token 写入仓库、命令历史或文档。

仓库根目录的 `.pypirc` 是本地凭据文件，受 `.gitignore` 管理，只用于本机发布，不能提交
到 GitHub。发布后仍应将其保存在安全位置。

## 正式 PyPI 发布

确认 TestPyPI 安装和导入成功后，再发布到 PyPI：

```bash
python -m twine upload dist/*
python -m pip install --upgrade stiffgwpy_fast
```

每次正式发布前都要递增 `pyproject.toml` 中的 `project.version`，并同步更新
`CHANGELOG.md`、`README.md` 和 `README_zh.md`。项目元数据中的联系邮箱为
`2966684515@qq.com`；PyPI 上传使用 API token，不使用邮箱密码。
