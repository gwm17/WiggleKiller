// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME srcdIsps_dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "include/DataStructs.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *DPPChannel_Dictionary();
   static void DPPChannel_TClassManip(TClass*);
   static void *new_DPPChannel(void *p = 0);
   static void *newArray_DPPChannel(Long_t size, void *p);
   static void delete_DPPChannel(void *p);
   static void deleteArray_DPPChannel(void *p);
   static void destruct_DPPChannel(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DPPChannel*)
   {
      ::DPPChannel *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::DPPChannel));
      static ::ROOT::TGenericClassInfo 
         instance("DPPChannel", "include/DataStructs.h", 11,
                  typeid(::DPPChannel), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &DPPChannel_Dictionary, isa_proxy, 4,
                  sizeof(::DPPChannel) );
      instance.SetNew(&new_DPPChannel);
      instance.SetNewArray(&newArray_DPPChannel);
      instance.SetDelete(&delete_DPPChannel);
      instance.SetDeleteArray(&deleteArray_DPPChannel);
      instance.SetDestructor(&destruct_DPPChannel);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DPPChannel*)
   {
      return GenerateInitInstanceLocal((::DPPChannel*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::DPPChannel*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *DPPChannel_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::DPPChannel*)0x0)->GetClass();
      DPPChannel_TClassManip(theClass);
   return theClass;
   }

   static void DPPChannel_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *DetectorHit_Dictionary();
   static void DetectorHit_TClassManip(TClass*);
   static void *new_DetectorHit(void *p = 0);
   static void *newArray_DetectorHit(Long_t size, void *p);
   static void delete_DetectorHit(void *p);
   static void deleteArray_DetectorHit(void *p);
   static void destruct_DetectorHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::DetectorHit*)
   {
      ::DetectorHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::DetectorHit));
      static ::ROOT::TGenericClassInfo 
         instance("DetectorHit", "include/DataStructs.h", 17,
                  typeid(::DetectorHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &DetectorHit_Dictionary, isa_proxy, 4,
                  sizeof(::DetectorHit) );
      instance.SetNew(&new_DetectorHit);
      instance.SetNewArray(&newArray_DetectorHit);
      instance.SetDelete(&delete_DetectorHit);
      instance.SetDeleteArray(&deleteArray_DetectorHit);
      instance.SetDestructor(&destruct_DetectorHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::DetectorHit*)
   {
      return GenerateInitInstanceLocal((::DetectorHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::DetectorHit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *DetectorHit_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::DetectorHit*)0x0)->GetClass();
      DetectorHit_TClassManip(theClass);
   return theClass;
   }

   static void DetectorHit_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *SabreDetector_Dictionary();
   static void SabreDetector_TClassManip(TClass*);
   static void *new_SabreDetector(void *p = 0);
   static void *newArray_SabreDetector(Long_t size, void *p);
   static void delete_SabreDetector(void *p);
   static void deleteArray_SabreDetector(void *p);
   static void destruct_SabreDetector(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SabreDetector*)
   {
      ::SabreDetector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::SabreDetector));
      static ::ROOT::TGenericClassInfo 
         instance("SabreDetector", "include/DataStructs.h", 22,
                  typeid(::SabreDetector), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &SabreDetector_Dictionary, isa_proxy, 4,
                  sizeof(::SabreDetector) );
      instance.SetNew(&new_SabreDetector);
      instance.SetNewArray(&newArray_SabreDetector);
      instance.SetDelete(&delete_SabreDetector);
      instance.SetDeleteArray(&deleteArray_SabreDetector);
      instance.SetDestructor(&destruct_SabreDetector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SabreDetector*)
   {
      return GenerateInitInstanceLocal((::SabreDetector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SabreDetector*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *SabreDetector_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::SabreDetector*)0x0)->GetClass();
      SabreDetector_TClassManip(theClass);
   return theClass;
   }

   static void SabreDetector_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *FPDetector_Dictionary();
   static void FPDetector_TClassManip(TClass*);
   static void *new_FPDetector(void *p = 0);
   static void *newArray_FPDetector(Long_t size, void *p);
   static void delete_FPDetector(void *p);
   static void deleteArray_FPDetector(void *p);
   static void destruct_FPDetector(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FPDetector*)
   {
      ::FPDetector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::FPDetector));
      static ::ROOT::TGenericClassInfo 
         instance("FPDetector", "include/DataStructs.h", 27,
                  typeid(::FPDetector), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &FPDetector_Dictionary, isa_proxy, 4,
                  sizeof(::FPDetector) );
      instance.SetNew(&new_FPDetector);
      instance.SetNewArray(&newArray_FPDetector);
      instance.SetDelete(&delete_FPDetector);
      instance.SetDeleteArray(&deleteArray_FPDetector);
      instance.SetDestructor(&destruct_FPDetector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FPDetector*)
   {
      return GenerateInitInstanceLocal((::FPDetector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::FPDetector*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *FPDetector_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::FPDetector*)0x0)->GetClass();
      FPDetector_TClassManip(theClass);
   return theClass;
   }

   static void FPDetector_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *CoincEvent_Dictionary();
   static void CoincEvent_TClassManip(TClass*);
   static void *new_CoincEvent(void *p = 0);
   static void *newArray_CoincEvent(Long_t size, void *p);
   static void delete_CoincEvent(void *p);
   static void deleteArray_CoincEvent(void *p);
   static void destruct_CoincEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::CoincEvent*)
   {
      ::CoincEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::CoincEvent));
      static ::ROOT::TGenericClassInfo 
         instance("CoincEvent", "include/DataStructs.h", 33,
                  typeid(::CoincEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &CoincEvent_Dictionary, isa_proxy, 4,
                  sizeof(::CoincEvent) );
      instance.SetNew(&new_CoincEvent);
      instance.SetNewArray(&newArray_CoincEvent);
      instance.SetDelete(&delete_CoincEvent);
      instance.SetDeleteArray(&deleteArray_CoincEvent);
      instance.SetDestructor(&destruct_CoincEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::CoincEvent*)
   {
      return GenerateInitInstanceLocal((::CoincEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::CoincEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *CoincEvent_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::CoincEvent*)0x0)->GetClass();
      CoincEvent_TClassManip(theClass);
   return theClass;
   }

   static void CoincEvent_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *ProcessedEvent_Dictionary();
   static void ProcessedEvent_TClassManip(TClass*);
   static void *new_ProcessedEvent(void *p = 0);
   static void *newArray_ProcessedEvent(Long_t size, void *p);
   static void delete_ProcessedEvent(void *p);
   static void deleteArray_ProcessedEvent(void *p);
   static void destruct_ProcessedEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ProcessedEvent*)
   {
      ::ProcessedEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::ProcessedEvent));
      static ::ROOT::TGenericClassInfo 
         instance("ProcessedEvent", "include/DataStructs.h", 38,
                  typeid(::ProcessedEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &ProcessedEvent_Dictionary, isa_proxy, 4,
                  sizeof(::ProcessedEvent) );
      instance.SetNew(&new_ProcessedEvent);
      instance.SetNewArray(&newArray_ProcessedEvent);
      instance.SetDelete(&delete_ProcessedEvent);
      instance.SetDeleteArray(&deleteArray_ProcessedEvent);
      instance.SetDestructor(&destruct_ProcessedEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ProcessedEvent*)
   {
      return GenerateInitInstanceLocal((::ProcessedEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ProcessedEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *ProcessedEvent_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::ProcessedEvent*)0x0)->GetClass();
      ProcessedEvent_TClassManip(theClass);
   return theClass;
   }

   static void ProcessedEvent_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_DPPChannel(void *p) {
      return  p ? new(p) ::DPPChannel : new ::DPPChannel;
   }
   static void *newArray_DPPChannel(Long_t nElements, void *p) {
      return p ? new(p) ::DPPChannel[nElements] : new ::DPPChannel[nElements];
   }
   // Wrapper around operator delete
   static void delete_DPPChannel(void *p) {
      delete ((::DPPChannel*)p);
   }
   static void deleteArray_DPPChannel(void *p) {
      delete [] ((::DPPChannel*)p);
   }
   static void destruct_DPPChannel(void *p) {
      typedef ::DPPChannel current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::DPPChannel

namespace ROOT {
   // Wrappers around operator new
   static void *new_DetectorHit(void *p) {
      return  p ? new(p) ::DetectorHit : new ::DetectorHit;
   }
   static void *newArray_DetectorHit(Long_t nElements, void *p) {
      return p ? new(p) ::DetectorHit[nElements] : new ::DetectorHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_DetectorHit(void *p) {
      delete ((::DetectorHit*)p);
   }
   static void deleteArray_DetectorHit(void *p) {
      delete [] ((::DetectorHit*)p);
   }
   static void destruct_DetectorHit(void *p) {
      typedef ::DetectorHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::DetectorHit

namespace ROOT {
   // Wrappers around operator new
   static void *new_SabreDetector(void *p) {
      return  p ? new(p) ::SabreDetector : new ::SabreDetector;
   }
   static void *newArray_SabreDetector(Long_t nElements, void *p) {
      return p ? new(p) ::SabreDetector[nElements] : new ::SabreDetector[nElements];
   }
   // Wrapper around operator delete
   static void delete_SabreDetector(void *p) {
      delete ((::SabreDetector*)p);
   }
   static void deleteArray_SabreDetector(void *p) {
      delete [] ((::SabreDetector*)p);
   }
   static void destruct_SabreDetector(void *p) {
      typedef ::SabreDetector current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SabreDetector

namespace ROOT {
   // Wrappers around operator new
   static void *new_FPDetector(void *p) {
      return  p ? new(p) ::FPDetector : new ::FPDetector;
   }
   static void *newArray_FPDetector(Long_t nElements, void *p) {
      return p ? new(p) ::FPDetector[nElements] : new ::FPDetector[nElements];
   }
   // Wrapper around operator delete
   static void delete_FPDetector(void *p) {
      delete ((::FPDetector*)p);
   }
   static void deleteArray_FPDetector(void *p) {
      delete [] ((::FPDetector*)p);
   }
   static void destruct_FPDetector(void *p) {
      typedef ::FPDetector current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::FPDetector

namespace ROOT {
   // Wrappers around operator new
   static void *new_CoincEvent(void *p) {
      return  p ? new(p) ::CoincEvent : new ::CoincEvent;
   }
   static void *newArray_CoincEvent(Long_t nElements, void *p) {
      return p ? new(p) ::CoincEvent[nElements] : new ::CoincEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_CoincEvent(void *p) {
      delete ((::CoincEvent*)p);
   }
   static void deleteArray_CoincEvent(void *p) {
      delete [] ((::CoincEvent*)p);
   }
   static void destruct_CoincEvent(void *p) {
      typedef ::CoincEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::CoincEvent

namespace ROOT {
   // Wrappers around operator new
   static void *new_ProcessedEvent(void *p) {
      return  p ? new(p) ::ProcessedEvent : new ::ProcessedEvent;
   }
   static void *newArray_ProcessedEvent(Long_t nElements, void *p) {
      return p ? new(p) ::ProcessedEvent[nElements] : new ::ProcessedEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_ProcessedEvent(void *p) {
      delete ((::ProcessedEvent*)p);
   }
   static void deleteArray_ProcessedEvent(void *p) {
      delete [] ((::ProcessedEvent*)p);
   }
   static void destruct_ProcessedEvent(void *p) {
      typedef ::ProcessedEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::ProcessedEvent

namespace ROOT {
   static TClass *vectorlEDetectorHitgR_Dictionary();
   static void vectorlEDetectorHitgR_TClassManip(TClass*);
   static void *new_vectorlEDetectorHitgR(void *p = 0);
   static void *newArray_vectorlEDetectorHitgR(Long_t size, void *p);
   static void delete_vectorlEDetectorHitgR(void *p);
   static void deleteArray_vectorlEDetectorHitgR(void *p);
   static void destruct_vectorlEDetectorHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<DetectorHit>*)
   {
      vector<DetectorHit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<DetectorHit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<DetectorHit>", -2, "vector", 214,
                  typeid(vector<DetectorHit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEDetectorHitgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<DetectorHit>) );
      instance.SetNew(&new_vectorlEDetectorHitgR);
      instance.SetNewArray(&newArray_vectorlEDetectorHitgR);
      instance.SetDelete(&delete_vectorlEDetectorHitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEDetectorHitgR);
      instance.SetDestructor(&destruct_vectorlEDetectorHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<DetectorHit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<DetectorHit>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEDetectorHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<DetectorHit>*)0x0)->GetClass();
      vectorlEDetectorHitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEDetectorHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEDetectorHitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DetectorHit> : new vector<DetectorHit>;
   }
   static void *newArray_vectorlEDetectorHitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<DetectorHit>[nElements] : new vector<DetectorHit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEDetectorHitgR(void *p) {
      delete ((vector<DetectorHit>*)p);
   }
   static void deleteArray_vectorlEDetectorHitgR(void *p) {
      delete [] ((vector<DetectorHit>*)p);
   }
   static void destruct_vectorlEDetectorHitgR(void *p) {
      typedef vector<DetectorHit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<DetectorHit>

namespace {
  void TriggerDictionaryInitialization_sps_dict_Impl() {
    static const char* headers[] = {
"include/DataStructs.h",
0
    };
    static const char* includePaths[] = {
"/usr/include/root",
"/home/gwm17/WiggleKiller/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "sps_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
struct __attribute__((annotate("$clingAutoload$include/DataStructs.h")))  DetectorHit;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
struct __attribute__((annotate("$clingAutoload$include/DataStructs.h")))  DPPChannel;
struct __attribute__((annotate("$clingAutoload$include/DataStructs.h")))  SabreDetector;
struct __attribute__((annotate("$clingAutoload$include/DataStructs.h")))  FPDetector;
struct __attribute__((annotate("$clingAutoload$include/DataStructs.h")))  CoincEvent;
struct __attribute__((annotate("$clingAutoload$include/DataStructs.h")))  ProcessedEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "sps_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "include/DataStructs.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"CoincEvent", payloadCode, "@",
"DPPChannel", payloadCode, "@",
"DetectorHit", payloadCode, "@",
"FPDetector", payloadCode, "@",
"ProcessedEvent", payloadCode, "@",
"SabreDetector", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("sps_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_sps_dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_sps_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_sps_dict() {
  TriggerDictionaryInitialization_sps_dict_Impl();
}
