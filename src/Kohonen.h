////////////////////////////////////////////////////////////////////////////
// Kohonen.h : main header file for KOHONEN.DLL
// See Kohonen.cxx for implementation.
// See Kohonen.ssc for implementation of the corresponding S function.
//
// KOHONEN_DLL_API assists exporting/importing functions defined in KOHONEN.DLL.
// It is for exporting if this header file is included in the source file (.c/.cxx) 
// that defines these functions, else importing.
//
// Symbols are exported if SP_CHAPTER_KOHONEN is defined, else imported.
// SP_CHAPTER_KOHONEN is automatically defined by the SPLUS Chapter Wizard. 
// Look at the setting..., C/C++ tap. 
////////////////////////////////////////////////////////////////////////////

#if !defined(S_KOHONEN_INCLUDED_)
#define S_KOHONEN_INCLUDED_

#undef KOHONEN_DLL_API
#if defined(WIN32) && defined(SP_CHAPTER_KOHONEN)
	#define KOHONEN_DLL_API(returnType) __declspec(dllexport) returnType __stdcall
#elif defined(WIN32)
	#define KOHONEN_DLL_API(returnType) __declspec(dllimport) returnType __stdcall
#else
	#define KOHONEN_DLL_API(returnType) extern returnType
#endif

#ifdef __cplusplus
extern "C" {
#endif

KOHONEN_DLL_API(void) get_coordinates(double*, long*, double*, double*, double*, long*);

KOHONEN_DLL_API(void) get_individual_ss(double*, long*, long*, double*, long*, long*, double*);

KOHONEN_DLL_API(void) get_inter_ss(double*, long*, long*, double*, double*);
KOHONEN_DLL_API(void) rows_dist_weights(double*, long*, long*, double*, long*, double*);
KOHONEN_DLL_API(void) KohonenC(double*, long*, long*, double*, long*, double*, long*, long*, double*, long*, double*, long*, double*, long*, double*, long*, double*);
KOHONEN_DLL_API(void) KohonenCBatch(double*, long*, long*, double*, long*, long*, long*, double*, double*, long*, long*, long*, double*, double*, double*);
KOHONEN_DLL_API(void) arrangedist(long*, long*, long*, long*, double*, double*, long*, long*, double*);

#ifdef __cplusplus
}
#endif

#endif //S_KOHONEN_INCLUDED_
