#pragma once

#include <algorithm>
#include <cctype>
#include <charconv>
#include <locale>
#include <string_view>
#include <vector>

namespace OpenABF::io_utils
{

/** @brief Compare two string_views, ignoring case */
static auto icase_compare(const std::string_view a, const std::string_view b)
    -> bool
{
    // not the same length
    if (a.length() != b.length()) {
        return false;
    }

    // iterate over the characters
    for (std::size_t i = 0; i < a.length(); ++i) {
        if (std::tolower(a[i]) != std::tolower(b[i])) {
            return false;
        }
    }

    // success
    return true;
}

/** @brief Left trim */
static auto trim_left(std::string_view s) -> std::string_view
{
    const auto& loc = std::locale();
    const auto* start = std::find_if_not(
        std::begin(s), std::end(s),
        [&loc](auto ch) -> bool { return std::isspace(ch, loc); });
    s.remove_prefix(std::distance(std::begin(s), start));
    return s;
}

/** @brief Right trim */
static auto trim_right(std::string_view s) -> std::string_view
{
    const auto& loc = std::locale();
    const auto* start =
        std::find_if_not(s.rbegin(), s.rend(), [&loc](auto ch) -> bool {
            return std::isspace(ch, loc);
        }).base();
    s.remove_suffix(std::distance(start, std::end(s)));
    return s;
}

/** @brief Trim from both ends */
static auto trim(std::string_view s) -> std::string_view
{
    s = trim_left(s);
    s = trim_right(s);
    return s;
}

/**
 * @brief Split a string by a delimiter
 *
 * When provided conflicting delimiters, the largest delimiter will take
 * precedence:
 *
 * ```{.cpp}
 * split("a->b->c", "-", "->");  // returns {"a", "b", "c"}
 * ```
 */
template <typename... Ds>
static auto split(std::string_view s, const Ds&... ds)
    -> std::vector<std::string_view>
{
    constexpr std::string_view DEFAULT_DELIM{" "};

    // Build delimiters list
    std::vector<std::string_view> delimiters;
    if (sizeof...(ds) > 0) {
        delimiters = {ds...};
    } else {
        delimiters.emplace_back(DEFAULT_DELIM);
    }

    // Get a list of all delimiter start pos and sizes
    std::vector<
        std::pair<std::string_view::size_type, std::string_view::size_type>>
        delimPos;
    for (const auto& delim : delimiters) {
        auto b = s.find(delim, 0);
        while (b != std::string_view::npos) {
            delimPos.emplace_back(b, delim.size());
            b = s.find(delim, b + delim.size());
        }
    }

    // Sort the delimiter start positions by first and largest
    std::sort(
        delimPos.begin(), delimPos.end(),
        [](const auto& l, const auto& r) { return l.second > r.second; });
    std::sort(
        delimPos.begin(), delimPos.end(),
        [](const auto& l, const auto& r) { return l.first < r.first; });

    // Split string
    std::vector<std::string_view> tokens;
    std::string_view::size_type begin{0};
    for (const auto& [end, size] : delimPos) {
        // ignore nested delimiters
        if (end < begin) {
            continue;
        }
        // get from begin to delim start
        if (auto t = s.substr(begin, end - begin); not t.empty()) {
            tokens.emplace_back(t);
        }
        begin = end + size;
    }
    if (auto t = s.substr(begin); not t.empty()) {
        tokens.emplace_back(t);
    }

    return tokens;
}

/**
 * @brief Convenience wrapper around std::to_chars for converting numerics to
 * std::string_view
 *
 * Useful during file writing operations when you're reusing a buffer, but
 * don't want to duplicate the error checking code of using std::to_chars.
 */
template <typename T>
auto to_string_view(const T& a, char* buf, const std::size_t& bufSize)
{
    auto res = std::to_chars(buf, buf + bufSize, a);
    if (res.ec != std::errc()) {
        throw std::runtime_error(std::make_error_code(res.ec).message());
    }
    return std::string_view(buf, res.ptr - buf);
}

/**
 * @brief Convert a string to a numeric type.
 *
 * A drop-in replacement for the `std:sto` family of functions which uses
 * `std::from_chars` for conversion. Like `std::sto`, throws exceptions when
 * conversion fails or if the converted value is out of range of the result
 * type.
 *
 * @throws std::invalid_argument If string cannot be converted to the result
 * type.
 * @throws std::result_out_of_range If converted value is out of range for the
 * result type.
 * @tparam T Requested numeric type
 * @tparam Args Parameter pack type
 * @param str Value to convert
 * @param args Extra parameters passed directly to `std::to_chars`
 * @return Converted value
 */
template <typename T, typename... Args>
auto to_numeric(const std::string_view str, Args... args) -> T
{
    T val;
    const auto* first = std::data(str);
    const auto* last = std::data(str) + std::size(str);
    auto [ptr, ec] = std::from_chars(first, last, val, args...);
    if (ec == std::errc::invalid_argument) {
        throw std::invalid_argument("Conversion could not be performed");
    }
    if (ec == std::errc::result_out_of_range) {
        throw std::out_of_range("Value out of range for the result type");
    }
    return val;
}

/**
 * @copybrief to_numeric
 *
 * Template specialization as fallback when the compiler does not support
 * `std::from_chars` for floating point types. Converts the input to a
 * `std::string` and passes to the appropriate `std::sto` function.
 */
template <>
inline auto to_numeric<float>(const std::string_view str) -> float
{
    return std::stof(std::string(str));
}

/** @copydoc to_numeric<float> */
template <>
inline auto to_numeric<double>(const std::string_view str) -> double
{
    return std::stod(std::string(str));
}

/** @copydoc to_numeric<float> */
template <>
inline auto to_numeric<long double>(const std::string_view str) -> long double
{
    return std::stold(std::string(str));
}
}  // namespace OpenABF::io_utils
