all:									libAliHLTTPCCAGPU.so libAliHLTTPCCAGPUOpenCL.so

clean:
										rm -f libAliHLTTPCCAGPU*.so AliHLTTPCCAGPUTracker*.o G__AliHLTTPCCAGPU*.o AliHLTTPCCAGPUTrackerNVCC.cu.cxx G__AliHLTTPCCAGPUAutoLinkDef*.h G__AliHLTTPCCAGPU*.h G__AliHLTTPCCAGPU*.cxx makefiles/opencl_compiler AliHLTTPCCAGPUTrackerOpenCLCode.*

libAliHLTTPCCAGPU.so:					AliHLTTPCCAGPUTrackerNVCC.o AliHLTTPCCAGPUTrackerBase.o G__AliHLTTPCCAGPUNVCC.o
										c++ -shared -L/usr/local/cuda/lib64 -L/opt/cuda/lib64 -L${ALICE_ROOT}/lib/tgt_${ALICE_TARGET} -L. -lcuda -lcudart -lAliHLTTPC -o $@ $^

libAliHLTTPCCAGPUOpenCL.so:				AliHLTTPCCAGPUTrackerOpenCL.o AliHLTTPCCAGPUTrackerBase.o G__AliHLTTPCCAGPUOpenCL.o AliHLTTPCCAGPUTrackerOpenCLCode.o
										c++ -shared -L$(AMDAPPSDKROOT)/lib/x86_64 -L${ALICE_ROOT}/lib/tgt_${ALICE_TARGET} -L. -lOpenCL -lAliHLTTPC -o $@ $^

AliHLTTPCCAGPUTrackerNVCC.o:			AliHLTTPCCAGPUTrackerNVCC.cu.cxx
										c++ -fPIC -DPACKAGE_TARNAME=\"alice-hlt\" -DPACKAGE_VERSION=\"35631\" -DPACKAGE_BUGREPORT=\"Matthias.Richter@ift.uib.no\" -DPACKAGE=\"alice-hlt\" -DVERSION=\"35631\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_DLFCN_H=1 -DLT_OBJDIR=\".libs/\" -DNDEBUG=1 -Duse_aliroot=1 -Duse_root=1 -DHAVE_HOMERREADER=1 -DHLT_SAMPLE=1 -DHLT_UTIL=1 -DHAVE_ALITPCRAWSTREAM_H=1 -DHLT_TPC=1 -DHAVE_NOT_TPCOFFLINE_REC=1 -DHAVE_TPC_MAPPING=1 -DHAVE_ALIALTRODECODER_H=1 -DHLT_RCU=1 -DHAVE_ALICALORAWSTREAM=1 -DHLT_CALO=1 -DHAVE_ALICALORAWSTREAM=1 -DHLT_PHOS=1 -DHLT_EMCAL=1 -DHLT_TRD=1 -DHLT_FMD=1 -DHAVE_ALIMPEXMAP_H=1 -DHAVE_ALIMUONTRIGGERIO_H=1 -DHLT_MUON=1 -DHLT_TRIGGER=1 -DHLT_GLOBAL=1 -DHLT_JET=1 -DHAVE_ALIITSCOMPRESSRAWDATASDD_H=1 -DHLT_ITS=1 -DHLT_COMP=1 -DMODULE=AliHLTTPC -W -Wall -Wshadow -Wno-effc++ -DROOTVERSION=\"5.25/02\" -DALIROOTVERSION=\"Unknown\" -O2 -DBUILD_GPU -c AliHLTTPCCAGPUTrackerNVCC.cu.cxx -o AliHLTTPCCAGPUTrackerNVCC.o

AliHLTTPCCAGPUTrackerBase.o:			AliHLTTPCCAGPUTrackerBase.cxx
										c++ -fPIC -DPACKAGE_TARNAME=\"alice-hlt\" -DPACKAGE_VERSION=\"35631\" -DPACKAGE_BUGREPORT=\"Matthias.Richter@ift.uib.no\" -DPACKAGE=\"alice-hlt\" -DVERSION=\"35631\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_DLFCN_H=1 -DLT_OBJDIR=\".libs/\" -DNDEBUG=1 -Duse_aliroot=1 -Duse_root=1 -DHAVE_HOMERREADER=1 -DHLT_SAMPLE=1 -DHLT_UTIL=1 -DHAVE_ALITPCRAWSTREAM_H=1 -DHLT_TPC=1 -DHAVE_NOT_TPCOFFLINE_REC=1 -DHAVE_TPC_MAPPING=1 -DHAVE_ALIALTRODECODER_H=1 -DHLT_RCU=1 -DHAVE_ALICALORAWSTREAM=1 -DHLT_CALO=1 -DHAVE_ALICALORAWSTREAM=1 -DHLT_PHOS=1 -DHLT_EMCAL=1 -DHLT_TRD=1 -DHLT_FMD=1 -DHAVE_ALIMPEXMAP_H=1 -DHAVE_ALIMUONTRIGGERIO_H=1 -DHLT_MUON=1 -DHLT_TRIGGER=1 -DHLT_GLOBAL=1 -DHLT_JET=1 -DHAVE_ALIITSCOMPRESSRAWDATASDD_H=1 -DHLT_ITS=1 -DHLT_COMP=1 -DMODULE=AliHLTTPC -W -Weffc++ -Wall -Wshadow -DROOTVERSION=\"5.25/02\" -DALIROOTVERSION=\"Unknown\" -O2 -DBUILD_GPU -I${ALICE_ROOT}/HLT/BASE -I${ALICE_ROOT}/HLT/TPCLib/tracking-ca -I${ROOTSYS}/include -c $< -o $@

AliHLTTPCCAGPUTrackerOpenCL.o:			AliHLTTPCCAGPUTrackerOpenCL.cxx
										c++ -fPIC -DPACKAGE_TARNAME=\"alice-hlt\" -DPACKAGE_VERSION=\"35631\" -DPACKAGE_BUGREPORT=\"Matthias.Richter@ift.uib.no\" -DPACKAGE=\"alice-hlt\" -DVERSION=\"35631\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_DLFCN_H=1 -DLT_OBJDIR=\".libs/\" -DNDEBUG=1 -Duse_aliroot=1 -Duse_root=1 -DHAVE_HOMERREADER=1 -DHLT_SAMPLE=1 -DHLT_UTIL=1 -DHAVE_ALITPCRAWSTREAM_H=1 -DHLT_TPC=1 -DHAVE_NOT_TPCOFFLINE_REC=1 -DHAVE_TPC_MAPPING=1 -DHAVE_ALIALTRODECODER_H=1 -DHLT_RCU=1 -DHAVE_ALICALORAWSTREAM=1 -DHLT_CALO=1 -DHAVE_ALICALORAWSTREAM=1 -DHLT_PHOS=1 -DHLT_EMCAL=1 -DHLT_TRD=1 -DHLT_FMD=1 -DHAVE_ALIMPEXMAP_H=1 -DHAVE_ALIMUONTRIGGERIO_H=1 -DHLT_MUON=1 -DHLT_TRIGGER=1 -DHLT_GLOBAL=1 -DHLT_JET=1 -DHAVE_ALIITSCOMPRESSRAWDATASDD_H=1 -DHLT_ITS=1 -DHLT_COMP=1 -DMODULE=AliHLTTPC -W -Weffc++ -Wall -Wshadow -DROOTVERSION=\"5.25/02\" -DALIROOTVERSION=\"Unknown\" -O2 -DBUILD_GPU -I$(AMDAPPSDKROOT)/include -I${ALICE_ROOT}/HLT/BASE -I${ALICE_ROOT}/HLT/TPCLib/tracking-ca -I${ROOTSYS}/include -Imakefiles -Wno-write-strings -c $< -o $@


G__AliHLTTPCCAGPU%.cxx:					AliHLTTPCCAGPUTracker%.h G__AliHLTTPCCAGPUAutoLinkDef%.h
										rootcint -f $@ -c -Duse_aliroot -Duse_root -DROWHOUGHPARAMS -Duse_reconstruction -Duse_newio -DROOTVERSION=\"unchecked\" -DALIROOTVERSION=\"unchecked\" -D__ROOT__ -DUSE_ALILOG -DLINUX -DNDEBUG -D_MODULE_=\"HLT\" -D`uname` -DDATE_SYS=`uname` -Dlong32='int' -Dlong64='long long' -DdatePointer='long' -I${ROOTSYS}/include -pthread -m64 -DWITHXML -DWITHXML -DUSE_ROOT -DWITHXML -I${ALICE_ROOT}/HLT/BASE -I${ALICE_ROOT}/HLT/BASE/util -I${ALICE_ROOT}/HLT -I${ALICE_ROOT}/HLT/TPCLib -I${ALICE_ROOT}/HLT/TPCLib/tracking-ca $^

G__AliHLTTPCCAGPUAutoLinkDef%.h:		AliHLTTPCCAGPUTracker%.h
										echo '//automatically generated ROOT DICT definition' > $@
										echo '//!!! DO NOT EDIT THIS FILE !!!' >> $@
										echo '#ifdef __CINT__' >> $@
										echo '#pragma link off all globals;' >> $@
										echo '#pragma link off all classes;' >> $@
										echo '#pragma link off all functions;' >> $@
										echo "#pragma link C++ class $<+;" | sed "s/\.h//" >> $@
										echo "#pragma link C++ class AliHLTTPCCAGPUTrackerBase+;" >> $@
										echo '#endif' >> $@

G__AliHLTTPCCAGPU%.o:					G__AliHLTTPCCAGPU%.cxx
										g++ -DcudaError_t=int -Duse_aliroot -Duse_root -DROWHOUGHPARAMS -Duse_reconstruction -Duse_newio -DROOTVERSION=\"unchecked\" -DALIROOTVERSION=\"unchecked\" -D__ROOT__ -DUSE_ALILOG -DLINUX -DNDEBUG -DBUILD_GPU -D_MODULE_=\"HLT\" -I${ALICE_ROOT}/HLT/TPCLib -I${ALICE_ROOT}/HLT/TPCLib/tracking-ca -I${ALICE_ROOT}/HLT/BASE -c $< -o $@ -O -g -W -Wall -Weffc++ -fPIC -pipe -fmessage-length=0 -Wno-long-long -ansi -Dlinux -D`uname` -DDATE_SYS=`uname` -Dlong32='int' -Dlong64='long long' -DdatePointer='long' -I${ROOTSYS}/include -pthread -m64 -D__PHOSUTIL__ -D__EMCALUTIL__

AliHLTTPCCAGPUTrackerNVCC.cu.cxx:		AliHLTTPCCAGPUTrackerNVCC.cu
										nvcc --cuda --use_fast_math --maxrregcount 64 -O4 -Xptxas -v -Xptxas -O4 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 --compiler-options "-DPACKAGE_TARNAME=\"alice-hlt\" -DPACKAGE_VERSION=\"35631\" -DPACKAGE_BUGREPORT=\"Matthias.Richter@ift.uib.no\" -DPACKAGE=\"alice-hlt\" -DVERSION=\"35631\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_DLFCN_H=1 -DLT_OBJDIR=\".libs/\" -DNDEBUG=1 -Duse_aliroot=1 -Duse_root=1 -DHAVE_HOMERREADER=1 -DHLT_SAMPLE=1 -DHLT_UTIL=1 -DHAVE_ALITPCRAWSTREAM_H=1 -DHLT_TPC=1 -DHAVE_NOT_TPCOFFLINE_REC=1 -DHAVE_TPC_MAPPING=1 -DHAVE_ALIALTRODECODER_H=1 -DHLT_RCU=1 -DHAVE_ALICALORAWSTREAM=1 -DHLT_CALO=1 -DHAVE_ALICALORAWSTREAM=1 -DHLT_PHOS=1 -DHLT_EMCAL=1 -DHLT_TRD=1 -DHLT_FMD=1 -DHAVE_ALIMPEXMAP_H=1 -DHAVE_ALIMUONTRIGGERIO_H=1 -DHLT_MUON=1 -DHLT_TRIGGER=1 -DHLT_GLOBAL=1 -DHLT_JET=1 -DHAVE_ALIITSCOMPRESSRAWDATASDD_H=1 -DHLT_ITS=1 -DHLT_COMP=1 -DMODULE=AliHLTTPC -IRCU -W -Wall -Wshadow -DROOTVERSION=\"5.25/02\" -DALIROOTVERSION=\"Unknown\" -O2 -DBUILD_GPU -I${ALICE_ROOT}/HLT/BASE -I${ALICE_ROOT}/HLT/TPCLib/tracking-ca -I${ROOTSYS}/include" -I. $< --output-file $@

AliHLTTPCCAGPUTrackerOpenCLCode.o:		AliHLTTPCCAGPUTrackerOpenCLCode.bin
										gcc -c makefiles/include.S -o $@

AliHLTTPCCAGPUTrackerOpenCLCode.bin:	AliHLTTPCCAGPUTrackerOpenCL.cl makefiles/opencl_compiler
										makefiles/opencl_compiler -output-file $@ AliHLTTPCCAGPUTrackerOpenCL.cl -- -I. -I${ALICE_ROOT}/HLT/BASE -I${ALICE_ROOT}/HLT/TPCLib/tracking-ca -I${ROOTSYS}/include -x clc++

makefiles/opencl_compiler:				makefiles/makefile_opencl_compiler.cpp
										c++ $< -o $@ -I$(AMDAPPSDKROOT)/include -L$(AMDAPPSDKROOT)/lib/x86_64 -lOpenCL
