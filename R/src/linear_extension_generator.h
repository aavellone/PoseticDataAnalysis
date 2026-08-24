/**
 * @file linear_extension_generator.h
 * @brief Header file for the abstract base class LinearExtensionGenerator.
 */

#pragma once

#include <cstdint>
#include <string>

// Assuming these headers exist in your project structure
#include "linear_extension.h"
#include "my_exception.h"

/**
 * @class LinearExtensionGenerator
 * @brief Abstract base class for generating linear extensions.
 *
 * @details This class provides the fundamental interface and protected state
 * for all linear extension generators in the HPC environment. It manages
 * the core buffer and generation tracking.
 */
class LinearExtensionGenerator {
protected:
    LinearExtension current_linear_extension_; ///< Buffer for the current linear extension.
    std::uint64_t n_elements_;                 ///< Number of elements in the linear extension.
    std::uint64_t current_number_le_ = 0;      ///< Counter for the generated extensions.
    bool started_ = false;                     ///< Flag indicating if generation has started.
    
public:
    /**
     * @brief Constructs the base generator and allocates the extension buffer.
     *
     * @param size The number of elements the linear extension will contain.
     */
    explicit LinearExtensionGenerator(std::uint64_t size)
    : current_linear_extension_(size),
    n_elements_(size) {}
    
    /**
     * @brief Virtual default destructor.
     */
    virtual ~LinearExtensionGenerator() = default;
    
    /**
     * @brief Initializes the generator.
     *
     * @param id Optional identifier or starting parameter for the generation.
     */
    virtual void Start(std::uint64_t id = 0) = 0;
    
    /**
     * @brief Advances to the next linear extension.
     */
    virtual void Next() = 0;
    
    /**
     * @brief Checks if there is a subsequent linear extension available.
     *
     * @return true if another extension can be generated, false otherwise.
     */
    [[nodiscard]] virtual bool HasNext() = 0;
    
    /**
     * @brief Returns the total number of linear extensions that can be generated.
     *
     * @return The theoretical maximum number of generated extensions.
     */
    [[nodiscard]] virtual std::uint64_t NumberOfLe() const = 0;
    
    /**
     * @brief Gets the number of linear extensions generated so far.
     *
     * @return The current value of the internal counter.
     */
    [[nodiscard]] std::uint64_t CurrentNumberOfLe() const noexcept {
        return current_number_le_;
    }
    
    /**
     * @brief Retrieves the current linear extension.
     *
     * @return A constant reference to the underlying LinearExtension object.
     * @throws MyException If the generator has not been started yet.
     */
    [[nodiscard]] const LinearExtension& Get() const {
        if (!started_) {
            throw MyException("LEG error: not started yet!");
        }
        return current_linear_extension_;
    }
    
    /**
     * @brief Gets the size (number of elements) of the linear extension.
     *
     * @return The allocated size of the current extension.
     */
    [[nodiscard]] std::uint64_t LeSize() const noexcept {
        return current_linear_extension_.size();
    }
};

