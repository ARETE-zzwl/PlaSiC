from __future__ import annotations

from dataclasses import dataclass

# importlib.resources.files 用于访问 Python package 内部打包的资源文件。
#
# 这里不是直接用普通文件系统路径，例如：
#     open("resources/parameters.json")
#
# 而是从 Python package：
#     plasic_app.resources
#
# 中读取资源。
#
# 这样做通常更适合已经安装成 package / wheel 的 Python 项目，
# 因为资源文件不一定应该依赖当前工作目录（cwd）。
from importlib.resources import files

import json

# Any 用于表示字典中的值可以是任意 Python 类型。
from typing import Any


def _resource_json(name: str) -> dict[str, Any]:
    """
    从 plasic_app.resources package 中读取指定 JSON 资源文件，
    并将其解析成 Python 字典。

    参数
    ----------
    name:
        package 内部的资源文件名，例如：
        - "parameters.json"
        - "variables.json"

    返回
    ----------
    dict[str, Any]
        JSON 文件解析后的顶层字典。

    说明
    ----------
    函数名前面的下划线 `_resource_json` 表示：
    这是模块内部使用的辅助函数，通常不作为该模块的公共 API。
    """

    # files("plasic_app.resources")
    # 获取 plasic_app.resources 这个 Python package 对应的资源根目录。
    resource = files("plasic_app.resources").joinpath(name)

    return json.loads(resource.read_text(encoding="utf-8"))


@dataclass(frozen=True)
class ParameterCatalog:
    """
    参数目录（parameter catalog）。

    这个类本身主要承担两个职责：

    1. 保存 parameters.json 中加载出来的参数元数据；
    2. 提供一些方便查询这些参数的方法。
    """

    # parameters.json 的 schema 版本号。
    schema_version: int

    # 参数元数据集合。
    # 每个参数都以 dict[str, Any] 表示；
    # 整体保存为 tuple，而不是 list。
    parameters: tuple[dict[str, Any], ...]

    @classmethod
    def load(cls) -> "ParameterCatalog":
        """
        从内置 parameters.json 资源加载 ParameterCatalog。

        这是一个 classmethod，因此调用方式通常是：

            catalog = ParameterCatalog.load()

        而不需要先创建 ParameterCatalog 实例。
        """

        # 从：
        #     plasic_app.resources/parameters.json
        # 中读取完整 JSON 数据。
        payload = _resource_json("parameters.json")

        # 根据 JSON 内容构造当前类实例。
        return cls(
            schema_version=int(payload["schema_version"]),
            parameters=tuple(payload["parameters"]),
        )

    def exposed(self) -> tuple[dict[str, Any], ...]:
        """
        返回“对外暴露”的参数集合。

        某些参数虽然存在于总目录中，但并不应该直接暴露给外部 UI

        当前会排除以下三种 status：

        - catalogued_not_exposed
        - derived_not_exposed
        - toolbar_managed

        返回结果仍然是 tuple。
        """

        # 这里使用 generator expression 逐个遍历 self.parameters，
        # 然后由 tuple(...) 将筛选结果一次性转换成 tuple。
        return tuple(
            item
            for item in self.parameters
            if item.get("status")
            not in {
                "catalogued_not_exposed",
                "derived_not_exposed",
                "toolbar_managed",
            }
        )

    def by_name(self, name: str) -> dict[str, Any]:
        """
        根据参数名称查找单个参数定义。
        """

        # 按 parameters 中原始顺序逐个检查参数。
        for parameter in self.parameters:
            if parameter["name"] == name:
                return parameter

        raise KeyError(name)

    def defaults(self, *, exposed_only: bool = True) -> dict[str, Any]:
        """
        生成 parameter name -> default value 的映射。

        参数
        ----------
        exposed_only:
            是否只包含 exposed 参数。

            True（默认）：
                只处理 self.exposed() 返回的参数。

            False：
                处理 self.parameters 中的全部参数。

        返回示例
        ----------
        假设参数元数据类似：

            {
                "name": "n_days",
                "default": 365
            }

        则返回字典中会包含：

            {
                "n_days": 365
            }

        如果某个参数没有 "default" 字段，则其值会是 None。
        """
        source = self.exposed() if exposed_only else self.parameters

        return {item["name"]: item.get("default") for item in source}


@dataclass(frozen=True)
class VariableCatalog:
    """
    变量目录（variable catalog）。

    与 ParameterCatalog 类似，它用于保存 variables.json 中定义的
    variable 元数据，并提供基于 id / name 的查询接口。

    VariableCatalog 比 ParameterCatalog 多了一个 protocol_version 字段，
    表明 variable 定义除了自身 schema 版本外，
    还关联到某个 protocol 版本。
    """

    # variables.json 自身的数据结构/schema 版本。
    schema_version: int

    # variables 所对应的协议版本。
    protocol_version: int

    # 所有 variable 的元数据集合。
    variables: tuple[dict[str, Any], ...]

    @classmethod
    def load(cls) -> "VariableCatalog":
        """
        从内置 variables.json 资源文件加载 VariableCatalog。

        调用示例：

            catalog = VariableCatalog.load()
        """

        payload = _resource_json("variables.json")

        return cls(
            schema_version=int(payload["schema_version"]),
            protocol_version=int(payload["protocol_version"]),
            variables=tuple(payload["variables"]),
        )

    def by_id(self, variable_id: int) -> dict[str, Any]:
        """
        根据 variable id 查找变量定义。
        """

        for variable in self.variables:
            if int(variable["id"]) == variable_id:
                return variable

        raise KeyError(variable_id)

    def by_name(self, name: str) -> dict[str, Any]:

        for variable in self.variables:

            if variable["name"] == name:
                return variable

        raise KeyError(name)
