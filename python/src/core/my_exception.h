#pragma once

#include <format>
#include <source_location>
#include <stdexcept>
#include <string>

// ================================================================
// MyException
// ================================================================

class MyException : public std::runtime_error {
public:
    MyException(const std::string& message,
                const std::source_location loc = std::source_location::current())
    : std::runtime_error(
                         std::format("{}:{}: {}", loc.file_name(), loc.line(), message)) {}
};

// ================================================================
// MyWarning
// ================================================================

class MyWarning : public std::runtime_error {
public:
    MyWarning(const std::string& message,
              const std::source_location loc = std::source_location::current())
    : std::runtime_error(
                         std::format("{}:{}: {}", loc.file_name(), loc.line(), message)) {}
};
