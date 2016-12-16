#pragma once

#include <array>

namespace ibmisc {

// https://gist.github.com/klmr/2775736
// Might be in upcoming C++ standard
// http://en.cppreference.com/w/cpp/experimental/make_array
template <typename... T>
constexpr auto make_array(T&&... values) ->
        std::array<
            typename std::decay<
                typename std::common_type<T...>::type>::type,
            sizeof...(T)> {
    return std::array<
        typename std::decay<
            typename std::common_type<T...>::type>::type,
        sizeof...(T)>{std::forward<T>(values)...};
}

}
