/** Copyright (c) 2023 Jim-shop
 * bp is licensed under Mulan PubL v2.
 * You can use this software according to the terms and conditions of the Mulan PubL v2.
 * You may obtain a copy of Mulan PubL v2 at:
 *          http://license.coscl.org.cn/MulanPubL-2.0
 * THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND,
 * EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT,
 * MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
 * See the Mulan PubL v2 for more details.
 */
#include <string>
#include <random>
#include <fstream>
#include <iostream>
#include <chrono>
#include <random>
#include <functional>

#include "nlohmann/json.hpp"

#include "bp_decoder/bp_decoder.hpp"
#include "sparse_matrix/sparse_matrix.hpp"

using ::std::operator""s;
using ::std::operator""sv;

template <typename T>
::std::ostream &operator<<(::std::ostream &stream, ::std::vector<T> const &vec)
{
    stream << '[';
    if (!vec.empty())
    {
        auto it{vec.begin()};
        stream << *it;
        for (; it != vec.end(); it++)
            stream << ", " << *it;
    }
    stream << ']';
    return stream;
}

template <>
::std::ostream &operator<<(::std::ostream &stream, ::std::vector<uint8_t> const &vec)
{
    stream << '[';
    if (!vec.empty())
    {
        auto it{vec.begin()};
        stream << (int)*it;
        for (; it != vec.end(); it++)
            stream << ", " << (int)*it;
    }
    stream << ']';
    return stream;
}

class Config
{
public: // data
    int random_seed;
    int target_runs;
    ::bp_decoder::BpDecoder::Method bp_method;
    double bit_error_rate;
    int max_iter;
    ::std::string hx_alist;

public: // apis
    auto &from_json(::std::string const &config_file)
    {
        try
        {
            ::std::ifstream stream{config_file};
            ::nlohmann::json json;
            stream >> json;

            auto input_seed = json.at("random_seed"sv).get<int>();
            this->random_seed = input_seed < 0 ? ::std::random_device{}()
                                               : this->random_seed = input_seed;

            auto input_targetruns = json.at("target_runs"sv).get<int>();
            this->target_runs = input_targetruns;

            auto input_bpmethod = json.at("bp_method"sv).get<::std::string>();
            this->bp_method = input_bpmethod == "min_sum"s ? bp_decoder::BpDecoder::Method::MIN_SUM
                                                           : bp_decoder::BpDecoder::Method::PRODUCT_SUM;

            auto input_biterrorrate = json.at("bit_error_rate"sv).get<double>();
            this->bit_error_rate = input_biterrorrate;

            auto input_maxiter = json.at("max_iter"sv).get<int>();
            this->max_iter = input_maxiter;

            auto hx_alist = json.at("hx_alist"sv).get<::std::string>();
            this->hx_alist = hx_alist;
        }
        catch (::nlohmann::json::parse_error const &err)
        {
            ::std::cerr << "语法错误或指定JSON文件名无法读取\n\n"sv
                        << err.what();
            ::std::exit(1);
        }
        catch (::nlohmann::json::out_of_range const &err)
        {
            ::std::cerr << "JSON文件未包含所需字段\n\n"sv
                        << err.what();
            ::std::exit(1);
        }

        return *this;
    }
    auto to_json() const
    {
        return ::nlohmann::json{
            {"random_seed"sv, this->random_seed},
            {"target_runs"sv, this->target_runs},
            {"bp_method"sv, this->bp_method == bp_decoder::BpDecoder::Method::MIN_SUM ? "min_sum"sv : "product_sum"sv},
            {"bit_error_rate"sv, this->bit_error_rate},
            {"max_iter"sv, this->max_iter}};
    }
};

class RandBitGen
{
    ::std::mt19937 rand_example;
    uint32_t threshold;

public:
    RandBitGen(uint32_t random_seed, double bit_error_rate)
        : rand_example{random_seed},
          threshold{static_cast<uint32_t>(bit_error_rate * UINT32_MAX)} {}
    uint8_t operator()() { return rand_example() < threshold; }
};

class Test
{
private: // vars
    ::RandBitGen randBitGen;
    ::bp_decoder::BpDecoder bpDecoder;
    int target_runs;
    ::sparse_matrix::Mod2SparseMatrix hx;

    ::std::vector<uint8_t> bit_error;

private: // utils
    // 返回是否生成了错误
    bool generateError(::std::vector<uint8_t> &arr)
    {
        bool has_error = false;
        for (auto &item : arr)
        {
            item = randBitGen();
            has_error |= item;
        }
        return has_error;
    }

public: // apis
    Test(::Config const &config)
        : randBitGen{static_cast<uint32_t>(config.random_seed), config.bit_error_rate},
          bpDecoder{config.bp_method, config.bit_error_rate, config.max_iter},
          target_runs{config.target_runs}
    {
        ::std::ifstream{config.hx_alist} >> this->hx;
        // ::std::cout << this->hx;
        this->bit_error.resize(hx.col);
        this->run();
    }
    void run()
    {
        auto bp_success_count{0ULL};
        for (auto curr{0}; curr < target_runs; curr++)
        {
            auto has_error = this->generateError(bit_error);
            if (!has_error)
                bp_success_count++; // treated as succeeded
            else
            {
                auto [run_iter, converge, log_prob_ratios, bp_decoding] = bpDecoder.run(this->hx, this->bit_error);
                // @todo
            }
        }
    }
};

inline static auto parseCommandLine(int argc, char *argv[])
{
    if (argc != 2)
    {
        ::std::cerr << "Json file input should be exactly one. \n"sv;
        exit(1);
    }
    return Config().from_json(argv[1]);
}

inline static auto timeit(::std::function<void()> f)
{
    auto start{::std::chrono::steady_clock::now()};
    f();
    auto end{::std::chrono::steady_clock::now()};
    return ::std::chrono::duration<double>(end - start);
}

int main(int argc, char *argv[])
{
    auto config = parseCommandLine(argc, argv);

    ::Test test(config);
    auto duration = timeit([&test]()
                           { test.run(); });
    ::std::cout << "Sim running time: "sv << duration;

    return 0;
}