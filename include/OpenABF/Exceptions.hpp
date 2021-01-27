#pragma once

#include <exception>

namespace OpenABF
{

/** @brief Solver exception */
class SolverException : std::exception
{
public:
    /** @brief Constructor with message */
    explicit SolverException(const char* msg) : msg_{msg} {}
    /** @brief Constructor with message */
    explicit SolverException(std::string msg) : msg_{std::move(msg)} {}
    /** @brief Get the exception message */
    const char* what() const noexcept override { return msg_.c_str(); }
private:
    /** Exception message */
    std::string msg_;
};

}  // namespace OpenABF