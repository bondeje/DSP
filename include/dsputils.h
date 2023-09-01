#ifdef DLL_EXPORT
#ifdef __linux__
#define EXPORT __attribute__((visibility("default")))
// else on Windows...use dllexport
#elif WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT 
#endif
#else
#define EXPORT 
#endif // DLL_EXPORT