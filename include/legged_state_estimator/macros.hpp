#ifndef INEKF_MACROS_HPP_
#define INEKF_MACROS_HPP_


#define INEKF_USE_DEFAULT_DESTTUCTOR(ClassName) ~ClassName() = default;
#define INEKF_USE_DEFAULT_COPY_CONSTRUCTOR(ClassName) ClassName(const ClassName&) = default;
#define INEKF_USE_DEFAULT_COPY_ASSIGN_OPERATOR(ClassName) ClassName& operator=(const ClassName&) = default;
#define INEKF_USE_DEFAULT_MOVE_CONSTRUCTOR(ClassName) ClassName(ClassName&&) noexcept = default;
#define INEKF_USE_DEFAULT_MOVE_ASSIGN_OPERATOR(ClassName) ClassName& operator=(ClassName&&) noexcept = default;


#endif // INEKF_MACROS_HPP_