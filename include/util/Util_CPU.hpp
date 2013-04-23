#pragma once


// Enable performance options in RELEASE mode
#if defined (NDEBUG) || defined (TREE_FMM_NDEBUG)

#ifndef TREE_FMM_INLINE
#define TREE_FMM_INLINE inline
#endif

#ifndef TREE_FMM_CHECK
#define TREE_FMM_CHECK 0
#endif

// Disable performance options in DEBUG mode
#else

#ifndef TREE_FMM_INLINE
#define TREE_FMM_INLINE
#endif

#ifndef TREE_FMM_CHECK
#define TREE_FMM_CHECK 1
#endif

#endif
