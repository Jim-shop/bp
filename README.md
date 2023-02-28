# BP-Decoder

C++ 实现的信念传播译码器。


## 构建

推荐使用 CMake 进行构建。 
在项目根目录下：
```shell
cmake -S . -B build
cmake --build build
```

构建结束会在 `build/` 中得到两个构建产物：

- `sim.exe`

    是测试（仿真）程序，使用方法见下文。

- `bp.wheel`

    是 Python 模块，使用方法见下文。


### 运行（仿真）

示例：
```shell
cd build
./sim ../data/test.json
```
需要向仿真程序传入配置文件(JSON)的路径。

JSON 示例见 [data/test.json](data/test.json)，
语法如下：
```json
{
    "random_seed": <int>, // 若为负数则随机生成一个随机种子
    "target_runs": <int>, // 仿真运行的次数
    "bp_method": <str>, // [ "min_sum" | "product_sum" ]
    "bit_error_rate": <double>, // [0, 1] 之间的浮点数
    "max_iter": <int>, // BP最大迭代次数
    "hx_alist": "../data/test.alist", // 输入校验矩阵
    "output_path": "../data/output/" // 输出路径
}
```

### Python 模块

@TODO

## 项目依赖

- [nlohmann/json](https://github.com/nlohmann/json)
- [alist format](http://www.inference.org.uk/mackay/codes/alist.html)

## 许可证

[木兰公共许可证，第二版](https://license.coscl.org.cn/MulanPubL-2.0)
