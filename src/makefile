INLINE_LIMIT=-finline-limit=100000
#INLINE_WARN=-Winline
STATIC=-static

######### Profile
#CLFLAGS=-pg -g $(STATIC) -DNDEBUG $(INLINE_WARN) $(INLINE_LIMIT) 

######### Profile & optimize
#CLFLAGS=-pg -g -O3 $(STATIC) -DNDEBUG $(INLINE_WARN) $(INLINE_LIMIT) 

######## Optimize
#CLFLAGS=-O3 $(STATIC) $(INLINE_WARN) $(INLINE_LIMIT)

######## Optimize and turn off asserts
CLFLAGS=-O3 -DNDEBUG $(STATIC) $(INLINE_WARN) $(INLINE_LIMIT)

########## Debug
CLFLAGS=-g

BASEDIR=..

GPP=g++
FLEX=flex
BISON=bison

BIN=$(BASEDIR)/bin
OBJDIR=$(BASEDIR)/bin/obj
MF=makefile

SRCDIR=$(BASEDIR)/src
UTIL=$(SRCDIR)/util
PARSER=$(SRCDIR)/parser
LOGIC=$(SRCDIR)/logic
INFER=$(SRCDIR)/infer
LEARNWTS=$(SRCDIR)/learnwts
LEARNSTRUCT=$(SRCDIR)/learnstruct


#ALCHEMY2.0 LIFTED ALGORITHMS
INCLUDESUFFIX=/includes
SRCSUFFIX=/src
LIFTEDINFER=$(SRCDIR)/liftedinfer
LIFTEDINFERBRIDGE=$(LIFTEDINFER)/bridge
LIFTEDINFERCOMMON=$(LIFTEDINFER)/common
LIFTEDINFERBG=$(LIFTEDINFER)/liftedblockedgibbs
LIFTEDINFERIS=$(LIFTEDINFER)/liftedimportancesampling
LIFTEDINFERWM=$(LIFTEDINFER)/liftedweightedmodelcounting
LIFTEDINFERTEST=$(LIFTEDINFER)/tests

REPLACEFOLY=$(PARSER)/replacefoly.pl
REPLACEFOLCPP=$(PARSER)/replacefolcpp.pl

LIBCINC= -I. -I$(UTIL) -I$(LOGIC) -I$(PARSER) -I$(LEARNWTS) -I$(LEARNSTRUCT) \
	-I$(INFER) -I$(LIFTEDINFERBRIDGE)$(INCLUDESUFFIX) -I$(LIFTEDINFERCOMMON)$(INCLUDESUFFIX) \
	-I$(LIFTEDINFERBG)$(INCLUDESUFFIX) -I$(LIFTEDINFERIS)$(INCLUDESUFFIX) -I$(LIFTEDINFERWM)$(INCLUDESUFFIX) \
	-I$(LIFTEDINFERTEST)
LIBLINC= -L. -L$(UTIL) -L$(LOGIC) -L$(PARSER) -L$(LEARNWTS) -L$(LEARNSTRUCT) \
	-L$(INFER)
LIBLINK= -ldl

#CFLAGS=$(CLFLAGS) -Wno-deprecated -Wall $(LIBCINC)
CFLAGS=$(CLFLAGS) -w -Wno-deprecated -Wall $(LIBCINC)
LFLAGS=$(CLFLAGS) -Wno-deprecated -Wall $(LIBLINC)


#############################################################################
#add non-executables .cpp files to OBJS 
#if .cpp file reside in a dir other than the current one, add .o rule below
OBJS=$(OBJDIR)/arguments.o $(OBJDIR)/powerset.o \
	$(OBJDIR)/fol.o $(OBJDIR)/listobj.o \
	$(OBJDIR)/predicatetemplate.o $(OBJDIR)/functiontemplate.o \
	$(OBJDIR)/function.o $(OBJDIR)/predicate.o $(OBJDIR)/term.o \
	$(OBJDIR)/domain.o $(OBJDIR)/clause.o $(OBJDIR)/clausesampler.o \
	$(OBJDIR)/Polynomial.o $(OBJDIR)/random.o \
	$(OBJDIR)/logichelper.o \
	$(OBJDIR)/groundclause.o \
	$(OBJDIR)/inferutil.o \
	$(OBJDIR)/createhypercube.o \
	$(OBJDIR)/createhypercubebasic.o $(OBJDIR)/createhypercubedt.o \
	$(OBJDIR)/joinhypercube.o \
	$(OBJDIR)/refinehypercube.o \
	$(OBJDIR)/superpred.o $(OBJDIR)/superclause.o \
	$(OBJDIR)/node.o $(OBJDIR)/factor.o \
	$(OBJDIR)/liftedalgsconvertor.o \
	$(OBJDIR)/liftedalgshandler.o \
	$(OBJDIR)/decomposer.o \
	$(OBJDIR)/fileutils.o \
	$(OBJDIR)/heuristics.o \
	$(OBJDIR)/logdouble.o \
	$(OBJDIR)/lvrmln.o \
	$(OBJDIR)/mathutils.o \
	$(OBJDIR)/normalizer.o \
	$(OBJDIR)/propositionalresolution.o \
	$(OBJDIR)/queryupdater.o \
	$(OBJDIR)/samplingalgs.o \
	$(OBJDIR)/blockexactinference.o \
	$(OBJDIR)/clustercreator.o \
	$(OBJDIR)/lbgsampler.o \
	$(OBJDIR)/lvrcluster.o \
	$(OBJDIR)/ptpsampler.o \
	$(OBJDIR)/lisapproxinference.o \
	$(OBJDIR)/proposalconstructor.o \
	$(OBJDIR)/proposalstructure.o \
	$(OBJDIR)/cache.o \
	$(OBJDIR)/extensions.o \
	$(OBJDIR)/ptpsearch.o \
	$(OBJDIR)/resolver.o \
	$(OBJDIR)/sampler.o \
	$(OBJDIR)/splitter.o \
	$(OBJDIR)/unifier.o \
	$(OBJDIR)/parser.o \
	$(OBJDIR)/testmodelcount.o \
	$(OBJDIR)/wmconvertor.o \
	$(OBJDIR)/randomgenutil.o \
	$(OBJDIR)/testformulainfer.o \
	$(OBJDIR)/lvrptptreesampling.o \
	$(OBJDIR)/lvrptpsearchtree.o

	
#add .cpp files (including executables') to CPPFILES.
CPPFILES=$(UTIL)/arguments.cpp $(UTIL)/powerset.cpp $(UTIL)/random.cpp\
	$(PARSER)/fol.cpp $(PARSER)/listobj.cpp \
	$(LOGIC)/predicatetemplate.cpp $(LOGIC)/functiontemplate.cpp \
	$(LOGIC)/function.cpp $(LOGIC)/predicate.cpp $(LOGIC)/term.cpp \
	$(LOGIC)/domain.cpp  $(LOGIC)/clause.cpp $(LOGIC)/clausesampler.cpp\
	$(INFER)/groundclause.cpp \
	$(LEARNWTS)/learnwts.cpp $(INFER)/infer.cpp \
	$(LEARNSTRUCT)/learnstruct.cpp \
	$(LOGIC)/Polynomial.cpp\
	$(LOGIC)/logichelper.cpp\
	$(INFER)/inferutil.cpp \
	$(INFER)/createhypercube.cpp \
	$(INFER)/createhypercubebasic.cpp $(INFER)/createhypercubedt.cpp \
	$(INFER)/joinhypercube.cpp \
	$(INFER)/refinehypercube.cpp \
	$(INFER)/superpred.cpp $(INFER)/superclause.cpp \
	$(INFER)/node.cpp $(INFER)/factor.cpp \
	$(LIFTEDINFERBRIDGE)$(SRCSUFFIX)/liftedalgsconvertor.cpp \
	$(LIFTEDINFERBRIDGE)$(SRCSUFFIX)/liftedalgshandler.cpp \
	$(LIFTEDINFERCOMMON)$(SRCSUFFIX)/decomposer.cpp \
	$(LIFTEDINFERCOMMON)$(SRCSUFFIX)/fileutils.cpp \
	$(LIFTEDINFERCOMMON)$(SRCSUFFIX)/heuristics.cpp \
	$(LIFTEDINFERCOMMON)$(SRCSUFFIX)/logdouble.cpp \
	$(LIFTEDINFERCOMMON)$(SRCSUFFIX)/lvrmln.cpp \
	$(LIFTEDINFERCOMMON)$(SRCSUFFIX)/mathutils.cpp \
	$(LIFTEDINFERCOMMON)$(SRCSUFFIX)/normalizer.cpp \
	$(LIFTEDINFERCOMMON)$(SRCSUFFIX)/propositionalresolution.cpp \
	$(LIFTEDINFERCOMMON)$(SRCSUFFIX)/queryupdater.cpp \
	$(LIFTEDINFERCOMMON)$(SRCSUFFIX)/samplingalgs.cpp \
	$(LIFTEDINFERCOMMON)$(SRCSUFFIX)/randomgenutil.cpp \
	$(LIFTEDINFERCOMMON)$(SRCSUFFIX)/lvrptptreesampling.cpp \
	$(LIFTEDINFERCOMMON)$(SRCSUFFIX)/lvrptpsearchtree.cpp \
	$(LIFTEDINFERBG)$(SRCSUFFIX)/blockexactinference.cpp \
	$(LIFTEDINFERBG)$(SRCSUFFIX)/clustercreator.cpp \
	$(LIFTEDINFERBG)$(SRCSUFFIX)/lbgsampler.cpp \
	$(LIFTEDINFERBG)$(SRCSUFFIX)/lvrcluster.cpp \
	$(LIFTEDINFERBG)$(SRCSUFFIX)/ptpsampler.cpp \
	$(LIFTEDINFERIS)$(SRCSUFFIX)/lisapproxinference.cpp \
	$(LIFTEDINFERIS)$(SRCSUFFIX)/proposalconstructor.cpp \
	$(LIFTEDINFERIS)$(SRCSUFFIX)/proposalstructure.cpp \
	$(LIFTEDINFERWM)$(SRCSUFFIX)/cache.cpp \
	$(LIFTEDINFERWM)$(SRCSUFFIX)/extensions.cpp \
	$(LIFTEDINFERWM)$(SRCSUFFIX)/ptpsearch.cpp \
	$(LIFTEDINFERWM)$(SRCSUFFIX)/resolver.cpp \
	$(LIFTEDINFERWM)$(SRCSUFFIX)/sampler.cpp \
	$(LIFTEDINFERWM)$(SRCSUFFIX)/splitter.cpp \
	$(LIFTEDINFERWM)$(SRCSUFFIX)/unifier.cpp \
	$(LIFTEDINFERWM)$(SRCSUFFIX)/wmconvertor.cpp \
	$(LIFTEDINFERTEST)/parser.cpp \
	$(LIFTEDINFERTEST)/testmodelcount.cpp \
	$(LIFTEDINFERTEST)/testformulainfer.cpp
	
##############################################################################
#EXECUTABLES: add executables rules here

all: learnwts infer learnstruct liftedinfer runliftedinfertests

learnwts: $(OBJDIR)/learnwts.oo
	cp $(OBJDIR)/learnwts.oo $(BIN)/learnwts

infer: $(OBJDIR)/infer.oo
	cp $(OBJDIR)/infer.oo $(BIN)/infer

liftedinfer: $(OBJDIR)/liftedinfer.oo
	cp $(OBJDIR)/liftedinfer.oo $(BIN)/liftedinfer

learnstruct: $(OBJDIR)/learnstruct.oo
	cp $(OBJDIR)/learnstruct.oo $(BIN)/learnstruct

runliftedinfertests:$(OBJDIR)/runtestsmain.oo
	cp $(OBJDIR)/runtestsmain.oo $(BIN)/runliftedinfertests

$(OBJDIR)/learnwts.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/learnwts.o -c $(LEARNWTS)/learnwts.cpp

$(OBJDIR)/infer.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/infer.o -c $(INFER)/infer.cpp

$(OBJDIR)/liftedinfer.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/liftedinfer.o -c $(LIFTEDINFER)/liftedinfer.cpp
	
$(OBJDIR)/runtestsmain.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/runtestsmain.o -c $(LIFTEDINFERTEST)/runtestsmain.cpp

$(OBJDIR)/learnstruct.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/learnstruct.o -c $(LEARNSTRUCT)/learnstruct.cpp

%.oo: %.o $(OBJS) $(LIBS)
	echo $<
	$(GPP) $(LFLAGS) -o $@ $(OBJS) $< $(LIBLINK)


##############################################################################
#add add non-executables .o rules here

$(OBJDIR)/fol.o: $(PARSER)/fol.y $(PARSER)/follex.y
	$(REPLACEFOLY) $(PARSER)/fol.y $(BISON)
	$(FLEX) -o$(PARSER)/follex.cpp $(PARSER)/follex.y 
	$(BISON) $(PARSER)/fol.y;
	mv -f fol.tab.c  $(PARSER)/fol.cpp
	$(REPLACEFOLCPP) $(PARSER)/fol.cpp
	$(GPP) $(CFLAGS) -o $(OBJDIR)/fol.o -c $(PARSER)/fol.cpp


#.o rule
#$(OBJDIR)/file.o:
#	$(GPP) $(CFLAGS) -o $(OBJDIR)/file.o -c dir/file.cpp


$(OBJDIR)/arguments.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/arguments.o -c $(UTIL)/arguments.cpp

$(OBJDIR)/powerset.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/powerset.o -c $(UTIL)/powerset.cpp

$(OBJDIR)/random.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/random.o -c $(UTIL)/random.cpp

$(OBJDIR)/listobj.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/listobj.o -c $(PARSER)/listobj.cpp 

$(OBJDIR)/predicatetemplate.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/predicatetemplate.o -c $(LOGIC)/predicatetemplate.cpp 

$(OBJDIR)/functiontemplate.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/functiontemplate.o -c $(LOGIC)/functiontemplate.cpp 

$(OBJDIR)/function.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/function.o -c $(LOGIC)/function.cpp 

$(OBJDIR)/predicate.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/predicate.o -c $(LOGIC)/predicate.cpp 

$(OBJDIR)/term.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/term.o -c $(LOGIC)/term.cpp

$(OBJDIR)/domain.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/domain.o -c $(LOGIC)/domain.cpp  

$(OBJDIR)/clause.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/clause.o -c $(LOGIC)/clause.cpp  

$(OBJDIR)/clausesampler.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/clausesampler.o -c $(LOGIC)/clausesampler.cpp  

$(OBJDIR)/groundclause.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/groundclause.o -c $(INFER)/groundclause.cpp

$(OBJDIR)/Polynomial.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/Polynomial.o -c $(LOGIC)/Polynomial.cpp
	
$(OBJDIR)/logichelper.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/logichelper.o -c $(LOGIC)/logichelper.cpp

$(OBJDIR)/inferutil.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/inferutil.o -c $(INFER)/inferutil.cpp

$(OBJDIR)/createhypercube.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/createhypercube.o -c $(INFER)/createhypercube.cpp

$(OBJDIR)/createhypercubebasic.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/createhypercubebasic.o -c $(INFER)/createhypercubebasic.cpp

$(OBJDIR)/createhypercubedt.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/createhypercubedt.o -c $(INFER)/createhypercubedt.cpp

$(OBJDIR)/joinhypercube.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/joinhypercube.o -c $(INFER)/joinhypercube.cpp

$(OBJDIR)/refinehypercube.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/refinehypercube.o -c $(INFER)/refinehypercube.cpp

$(OBJDIR)/superpred.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/superpred.o -c $(INFER)/superpred.cpp

$(OBJDIR)/superclause.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/superclause.o -c $(INFER)/superclause.cpp

$(OBJDIR)/node.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/node.o -c $(INFER)/node.cpp

$(OBJDIR)/factor.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/factor.o -c $(INFER)/factor.cpp

$(OBJDIR)/liftedalgsconvertor.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/liftedalgsconvertor.o -c $(LIFTEDINFERBRIDGE)$(SRCSUFFIX)/liftedalgsconvertor.cpp
$(OBJDIR)/liftedalgshandler.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/liftedalgshandler.o -c $(LIFTEDINFERBRIDGE)$(SRCSUFFIX)/liftedalgshandler.cpp
$(OBJDIR)/decomposer.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/decomposer.o -c $(LIFTEDINFERCOMMON)$(SRCSUFFIX)/decomposer.cpp
$(OBJDIR)/fileutils.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/fileutils.o -c $(LIFTEDINFERCOMMON)$(SRCSUFFIX)/fileutils.cpp
$(OBJDIR)/heuristics.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/heuristics.o -c $(LIFTEDINFERCOMMON)$(SRCSUFFIX)/heuristics.cpp
$(OBJDIR)/logdouble.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/logdouble.o -c $(LIFTEDINFERCOMMON)$(SRCSUFFIX)/logdouble.cpp
$(OBJDIR)/lvrmln.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/lvrmln.o -c $(LIFTEDINFERCOMMON)$(SRCSUFFIX)/lvrmln.cpp
$(OBJDIR)/mathutils.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/mathutils.o -c $(LIFTEDINFERCOMMON)$(SRCSUFFIX)/mathutils.cpp
$(OBJDIR)/normalizer.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/normalizer.o -c $(LIFTEDINFERCOMMON)$(SRCSUFFIX)/normalizer.cpp
$(OBJDIR)/propositionalresolution.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/propositionalresolution.o -c $(LIFTEDINFERCOMMON)$(SRCSUFFIX)/propositionalresolution.cpp
$(OBJDIR)/queryupdater.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/queryupdater.o -c $(LIFTEDINFERCOMMON)$(SRCSUFFIX)/queryupdater.cpp
$(OBJDIR)/samplingalgs.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/samplingalgs.o -c $(LIFTEDINFERCOMMON)$(SRCSUFFIX)/samplingalgs.cpp
$(OBJDIR)/randomgenutil.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/randomgenutil.o -c $(LIFTEDINFERCOMMON)$(SRCSUFFIX)/randomgenutil.cpp
$(OBJDIR)/lvrptptreesampling.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/lvrptptreesampling.o -c $(LIFTEDINFERCOMMON)$(SRCSUFFIX)/lvrptptreesampling.cpp
$(OBJDIR)/lvrptpsearchtree.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/lvrptpsearchtree.o -c $(LIFTEDINFERCOMMON)$(SRCSUFFIX)/lvrptpsearchtree.cpp
	
$(OBJDIR)/blockexactinference.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/blockexactinference.o -c $(LIFTEDINFERBG)$(SRCSUFFIX)/blockexactinference.cpp
$(OBJDIR)/clustercreator.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/clustercreator.o -c $(LIFTEDINFERBG)$(SRCSUFFIX)/clustercreator.cpp
$(OBJDIR)/lbgsampler.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/lbgsampler.o -c $(LIFTEDINFERBG)$(SRCSUFFIX)/lbgsampler.cpp
$(OBJDIR)/lvrcluster.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/lvrcluster.o -c $(LIFTEDINFERBG)$(SRCSUFFIX)/lvrcluster.cpp
$(OBJDIR)/ptpsampler.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/ptpsampler.o -c $(LIFTEDINFERBG)$(SRCSUFFIX)/ptpsampler.cpp
$(OBJDIR)/lisapproxinference.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/lisapproxinference.o -c $(LIFTEDINFERIS)$(SRCSUFFIX)/lisapproxinference.cpp
$(OBJDIR)/proposalconstructor.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/proposalconstructor.o -c $(LIFTEDINFERIS)$(SRCSUFFIX)/proposalconstructor.cpp
$(OBJDIR)/proposalstructure.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/proposalstructure.o -c $(LIFTEDINFERIS)$(SRCSUFFIX)/proposalstructure.cpp
$(OBJDIR)/cache.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/cache.o -c $(LIFTEDINFERWM)$(SRCSUFFIX)/cache.cpp
$(OBJDIR)/extensions.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/extensions.o -c $(LIFTEDINFERWM)$(SRCSUFFIX)/extensions.cpp
$(OBJDIR)/ptpsearch.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/ptpsearch.o -c $(LIFTEDINFERWM)$(SRCSUFFIX)/ptpsearch.cpp
$(OBJDIR)/resolver.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/resolver.o -c $(LIFTEDINFERWM)$(SRCSUFFIX)/resolver.cpp
$(OBJDIR)/sampler.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/sampler.o -c $(LIFTEDINFERWM)$(SRCSUFFIX)/sampler.cpp
$(OBJDIR)/splitter.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/splitter.o -c $(LIFTEDINFERWM)$(SRCSUFFIX)/splitter.cpp
$(OBJDIR)/unifier.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/unifier.o -c $(LIFTEDINFERWM)$(SRCSUFFIX)/unifier.cpp
$(OBJDIR)/parser.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/parser.o -c $(LIFTEDINFERTEST)/parser.cpp
$(OBJDIR)/testmodelcount.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/testmodelcount.o -c $(LIFTEDINFERTEST)/testmodelcount.cpp
$(OBJDIR)/wmconvertor.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/wmconvertor.o -c $(LIFTEDINFERWM)$(SRCSUFFIX)/wmconvertor.cpp
$(OBJDIR)/testformulainfer.o:
	$(GPP) $(CFLAGS) -o $(OBJDIR)/testformulainfer.o -c $(LIFTEDINFERTEST)/testformulainfer.cpp

###############################################################################


clean:
	@rm -f $(OBJDIR)/*.o $(OBJDIR)/*.oo

depend:
	@echo 'updating the dependencies for:'
	@echo '    ' $(CPPFILES)
	@( \
	< $(MF) sed -n '1,/^###.*SUDDEN DEATH/p'; \
	    echo '#' ; \
	    echo '# dependencies generated on: ' `date` ; \
	    echo '#' ; \
	    for i in $(CPPFILES); do \
	        echo -n "$(OBJDIR)/" ; \
	        $(GPP) -MM $(CFLAGS) $$i ; \
	        echo; \
	    done \
	) > $(MF).new
	@mv $(MF) $(MF).last
	@mv $(MF).new $(MF)

############ Everything below here is subject to SUDDEN DEATH #####
#
# dependencies generated on:  Sun May 28 07:35:03 EDT 2017
#
../bin/obj/arguments.o: ../src/util/arguments.cpp ../src/util/arguments.h

../bin/obj/powerset.o: ../src/util/powerset.cpp ../src/util/powerset.h \
 ../src/util/array.h ../src/util/random.h

../bin/obj/random.o: ../src/util/random.cpp ../src/util/random.h

../bin/obj/fol.o: ../src/parser/fol.cpp ../src/parser/fol.h ../src/logic/mln.h \
 ../src/logic/clause.h ../src/logic/domain.h ../src/util/equalstr.h \
 ../src/util/dualmap.h ../src/util/array.h ../src/util/random.h \
 ../src/util/equalstr.h ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicatetemplate.h ../src/util/array.h \
 ../src/logic/functiontemplate.h ../src/logic/predicate.h \
 ../src/logic/term.h ../src/util/hash.h ../src/util/hasharray.h \
 ../src/logic/function.h ../src/util/hashstring.h ../src/util/hasharray.h \
 ../src/util/hashlist.h ../src/util/arraysaccessor.h \
 ../src/util/powerset.h ../src/util/multdarray.h ../src/util/hashint.h \
 ../src/util/hash.h ../src/logic/database.h ../src/logic/groundpreds.h \
 ../src/infer/groundpredicate.h ../src/logic/domain.h \
 ../src/logic/predicate.h ../src/infer/groundclause.h \
 ../src/logic/clausesampler.h ../src/logic/truefalsegroundingsstore.h \
 ../src/util/random.h ../src/logic/clausehelper.h \
 ../src/logic/mlnhelper.h ../src/logic/database.h \
 ../src/parser/follex.cpp ../src/parser/folhelper.h \
 ../src/parser/strfifolist.h ../src/logic/function.h \
 ../src/logic/clause.h ../src/logic/clausefactory.h \
 ../src/parser/listobj.h ../src/util/timer.h ../src/parser/internals.h

../bin/obj/listobj.o: ../src/parser/listobj.cpp ../src/parser/listobj.h \
 ../src/util/array.h ../src/util/random.h ../src/util/hashstring.h \
 ../src/util/hasharray.h ../src/logic/domain.h ../src/util/equalstr.h \
 ../src/util/dualmap.h ../src/util/array.h ../src/util/equalstr.h \
 ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicatetemplate.h ../src/logic/functiontemplate.h \
 ../src/logic/predicate.h ../src/logic/term.h ../src/util/hash.h \
 ../src/util/hasharray.h ../src/logic/function.h \
 ../src/util/arraysaccessor.h

../bin/obj/predicatetemplate.o: ../src/logic/predicatetemplate.cpp \
 ../src/logic/predicatetemplate.h ../src/util/array.h \
 ../src/util/random.h ../src/logic/domain.h ../src/util/equalstr.h \
 ../src/util/dualmap.h ../src/util/array.h ../src/util/equalstr.h \
 ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/functiontemplate.h ../src/logic/predicate.h \
 ../src/logic/term.h ../src/util/hash.h ../src/util/hasharray.h \
 ../src/logic/function.h

../bin/obj/functiontemplate.o: ../src/logic/functiontemplate.cpp \
 ../src/logic/functiontemplate.h ../src/logic/predicatetemplate.h \
 ../src/util/array.h ../src/util/random.h ../src/logic/domain.h \
 ../src/util/equalstr.h ../src/util/dualmap.h ../src/util/array.h \
 ../src/util/equalstr.h ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicate.h ../src/logic/term.h ../src/util/hash.h \
 ../src/util/hasharray.h ../src/logic/function.h

../bin/obj/function.o: ../src/logic/function.cpp ../src/logic/function.h \
 ../src/logic/functiontemplate.h ../src/logic/predicatetemplate.h \
 ../src/util/array.h ../src/util/random.h ../src/util/hash.h \
 ../src/util/array.h ../src/logic/term.h ../src/logic/domain.h \
 ../src/util/equalstr.h ../src/util/dualmap.h ../src/util/equalstr.h \
 ../src/util/constdualmap.h ../src/util/strint.h ../src/logic/predicate.h \
 ../src/util/hasharray.h

../bin/obj/predicate.o: ../src/logic/predicate.cpp ../src/logic/predicate.h \
 ../src/logic/predicatetemplate.h ../src/util/array.h \
 ../src/util/random.h ../src/logic/term.h ../src/util/hash.h \
 ../src/util/array.h ../src/util/hasharray.h ../src/logic/function.h \
 ../src/logic/functiontemplate.h ../src/logic/clause.h \
 ../src/logic/domain.h ../src/util/equalstr.h ../src/util/dualmap.h \
 ../src/util/equalstr.h ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/util/hashstring.h ../src/util/hasharray.h ../src/util/hashlist.h \
 ../src/util/arraysaccessor.h ../src/util/powerset.h \
 ../src/util/multdarray.h ../src/util/hashint.h ../src/util/hash.h \
 ../src/logic/database.h ../src/logic/groundpreds.h \
 ../src/infer/groundpredicate.h ../src/logic/domain.h \
 ../src/logic/predicate.h ../src/infer/groundclause.h \
 ../src/logic/clausesampler.h ../src/logic/truefalsegroundingsstore.h \
 ../src/util/random.h ../src/logic/clausehelper.h

../bin/obj/term.o: ../src/logic/term.cpp ../src/logic/term.h ../src/util/array.h \
 ../src/util/random.h ../src/logic/function.h \
 ../src/logic/functiontemplate.h ../src/logic/predicatetemplate.h \
 ../src/util/hash.h ../src/util/array.h ../src/logic/domain.h \
 ../src/util/equalstr.h ../src/util/dualmap.h ../src/util/equalstr.h \
 ../src/util/constdualmap.h ../src/util/strint.h ../src/logic/predicate.h \
 ../src/util/hasharray.h

../bin/obj/domain.o: ../src/logic/domain.cpp ../src/logic/clause.h \
 ../src/logic/domain.h ../src/util/equalstr.h ../src/util/dualmap.h \
 ../src/util/array.h ../src/util/random.h ../src/util/equalstr.h \
 ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicatetemplate.h ../src/util/array.h \
 ../src/logic/functiontemplate.h ../src/logic/predicate.h \
 ../src/logic/term.h ../src/util/hash.h ../src/util/hasharray.h \
 ../src/logic/function.h ../src/util/hashstring.h ../src/util/hasharray.h \
 ../src/util/hashlist.h ../src/util/arraysaccessor.h \
 ../src/util/powerset.h ../src/util/multdarray.h ../src/util/hashint.h \
 ../src/util/hash.h ../src/logic/database.h ../src/logic/groundpreds.h \
 ../src/infer/groundpredicate.h ../src/logic/domain.h \
 ../src/logic/predicate.h ../src/infer/groundclause.h \
 ../src/logic/clausesampler.h ../src/logic/truefalsegroundingsstore.h \
 ../src/util/random.h ../src/logic/clausehelper.h ../src/logic/mln.h \
 ../src/logic/mlnhelper.h

../bin/obj/clause.o: ../src/logic/clause.cpp ../src/logic/clause.h \
 ../src/logic/domain.h ../src/util/equalstr.h ../src/util/dualmap.h \
 ../src/util/array.h ../src/util/random.h ../src/util/equalstr.h \
 ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicatetemplate.h ../src/util/array.h \
 ../src/logic/functiontemplate.h ../src/logic/predicate.h \
 ../src/logic/term.h ../src/util/hash.h ../src/util/hasharray.h \
 ../src/logic/function.h ../src/util/hashstring.h ../src/util/hasharray.h \
 ../src/util/hashlist.h ../src/util/arraysaccessor.h \
 ../src/util/powerset.h ../src/util/multdarray.h ../src/util/hashint.h \
 ../src/util/hash.h ../src/logic/database.h ../src/logic/groundpreds.h \
 ../src/infer/groundpredicate.h ../src/logic/domain.h \
 ../src/logic/predicate.h ../src/infer/groundclause.h \
 ../src/logic/clausesampler.h ../src/logic/truefalsegroundingsstore.h \
 ../src/util/random.h ../src/logic/clausehelper.h ../src/infer/mrf.h \
 ../src/util/timer.h ../src/logic/mln.h ../src/logic/mlnhelper.h \
 ../src/infer/groundpredicate.h ../src/infer/variable.h \
 ../src/util/util.h ../src/infer/superclause.h ../src/infer/mrf.h \
 ../src/infer/variable.h ../src/infer/hypercube.h \
 ../src/infer/inferutil.h ../src/logic/clause.h \
 ../src/infer/hypercubeoperations.h ../src/infer/hypercubereverseindex.h \
 ../src/infer/hypercuberefinement.h

../bin/obj/clausesampler.o: ../src/logic/clausesampler.cpp \
 ../src/logic/clausesampler.h ../src/logic/database.h ../src/util/hash.h \
 ../src/util/array.h ../src/util/random.h ../src/logic/groundpreds.h \
 ../src/logic/predicate.h ../src/logic/predicatetemplate.h \
 ../src/util/array.h ../src/logic/term.h ../src/util/hasharray.h \
 ../src/logic/function.h ../src/logic/functiontemplate.h \
 ../src/logic/domain.h ../src/util/equalstr.h ../src/util/dualmap.h \
 ../src/util/equalstr.h ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/util/arraysaccessor.h ../src/infer/groundpredicate.h \
 ../src/logic/domain.h ../src/logic/predicate.h \
 ../src/infer/groundclause.h ../src/logic/truefalsegroundingsstore.h \
 ../src/util/random.h ../src/logic/clause.h ../src/util/hashstring.h \
 ../src/util/hasharray.h ../src/util/hashlist.h ../src/util/powerset.h \
 ../src/util/multdarray.h ../src/util/hashint.h ../src/util/hash.h \
 ../src/logic/clausehelper.h

../bin/obj/groundclause.o: ../src/infer/groundclause.cpp ../src/infer/groundclause.h \
 ../src/util/hasharray.h ../src/util/random.h ../src/util/array.h \
 ../src/util/hash.h ../src/util/array.h ../src/infer/groundpredicate.h \
 ../src/logic/domain.h ../src/util/equalstr.h ../src/util/dualmap.h \
 ../src/util/equalstr.h ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicatetemplate.h ../src/logic/functiontemplate.h \
 ../src/logic/predicate.h ../src/logic/term.h ../src/logic/function.h \
 ../src/logic/predicate.h ../src/logic/clause.h ../src/logic/domain.h \
 ../src/util/hashstring.h ../src/util/hasharray.h ../src/util/hashlist.h \
 ../src/util/arraysaccessor.h ../src/util/powerset.h \
 ../src/util/multdarray.h ../src/util/hashint.h ../src/util/hash.h \
 ../src/logic/database.h ../src/logic/groundpreds.h \
 ../src/infer/groundpredicate.h ../src/logic/clausesampler.h \
 ../src/logic/truefalsegroundingsstore.h ../src/util/random.h \
 ../src/logic/clausehelper.h ../src/logic/mln.h ../src/logic/clause.h \
 ../src/logic/mlnhelper.h

../bin/obj/learnwts.o: ../src/learnwts/learnwts.cpp ../src/util/arguments.h \
 ../src/infer/inferenceargs.h ../src/infer/inferutil.h \
 ../src/util/array.h ../src/util/random.h ../src/util/hashint.h \
 ../src/util/hasharray.h ../src/util/hash.h ../src/util/array.h \
 ../src/logic/clause.h ../src/logic/domain.h ../src/util/equalstr.h \
 ../src/util/dualmap.h ../src/util/equalstr.h ../src/util/constdualmap.h \
 ../src/util/strint.h ../src/logic/predicatetemplate.h \
 ../src/logic/functiontemplate.h ../src/logic/predicate.h \
 ../src/logic/term.h ../src/util/hash.h ../src/util/hasharray.h \
 ../src/logic/function.h ../src/util/hashstring.h ../src/util/hashlist.h \
 ../src/util/arraysaccessor.h ../src/util/powerset.h \
 ../src/util/multdarray.h ../src/logic/database.h \
 ../src/logic/groundpreds.h ../src/infer/groundpredicate.h \
 ../src/logic/domain.h ../src/logic/predicate.h \
 ../src/infer/groundclause.h ../src/logic/clausesampler.h \
 ../src/logic/truefalsegroundingsstore.h ../src/util/random.h \
 ../src/logic/clausehelper.h ../src/learnwts/lbfgsb.h ../src/util/timer.h \
 ../src/learnwts/pseudologlikelihood.h ../src/logic/mln.h \
 ../src/logic/clause.h ../src/logic/mlnhelper.h \
 ../src/learnwts/indextranslator.h ../src/logic/Polynomial.h \
 ../src/learnwts/discriminativelearner.h ../src/infer/infer.h \
 ../src/util/util.h ../src/parser/fol.h ../src/logic/database.h \
 ../src/infer/mrf.h ../src/infer/groundpredicate.h \
 ../src/learnwts/learnwts.h ../src/infer/inferenceargs.h \
 ../src/infer/maxwalksat.h ../src/infer/sat.h ../src/infer/inference.h \
 ../src/logic/variablestate.h ../src/infer/mrf.h \
 ../src/logic/hvariablestate.h ../src/logic/Polynomial.h \
 ../src/learnwts/lbfgsp.h ../src/logic/logsigmoid.h \
 ../src/logic/logichelper.h ../src/infer/maxwalksatparams.h \
 ../src/infer/samplesatparams.h ../src/infer/mcsat.h ../src/infer/mcmc.h \
 ../src/infer/mcmcparams.h ../src/infer/mcsatparams.h \
 ../src/infer/unitpropagation.h ../src/infer/gibbssampler.h \
 ../src/infer/gibbsparams.h ../src/infer/convergencetest.h \
 ../src/util/meanvariance.h ../src/infer/gelmanconvergencetest.h \
 ../src/infer/simulatedtempering.h \
 ../src/infer/simulatedtemperingparams.h ../src/infer/bp.h \
 ../src/infer/bpparams.h ../src/infer/twowaymessage.h \
 ../src/infer/superclause.h ../src/infer/variable.h \
 ../src/infer/hypercube.h ../src/infer/hypercubeoperations.h \
 ../src/infer/hypercubereverseindex.h ../src/infer/hypercuberefinement.h \
 ../src/infer/auxfactor.h ../src/infer/factor.h ../src/infer/link.h \
 ../src/infer/node.h ../src/infer/superpred.h \
 ../src/infer/clausecounter.h ../src/infer/factorgraph.h \
 ../src/logic/clausefactory.h ../src/infer/hmcsat.h \
 ../src/infer/hmaxwalksat.h ../src/logic/logichelper.h \
 ../src/infer/maxwalksat.h ../src/infer/bp.h ../src/learnwts/learnwts.h \
 ../src/infer/mcsat.h ../src/infer/gibbssampler.h \
 ../src/infer/simulatedtempering.h \
 ../src/liftedinfer/bridge/includes/lbginference.h ../src/infer/mcmc.h \
 ../src/infer/mcsatparams.h \
 ../src/liftedinfer/bridge/includes/liftedalgsconvertor.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/bridge/includes/liftedalgshandler.h \
 ../src/liftedinfer/common/includes/lvparams.h ../src/infer/mcmcparams.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/clustercreator.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/ptpsampler.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/blockexactinference.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/common/includes/heuristics.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/lvrnormpropagation.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/randomgenutil.h \
 ../src/liftedinfer/common/includes/lvrsingletonnormpropagation.h \
 ../src/liftedinfer/common/includes/lvrptptreesampling.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/lvrcluster.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/gibbsprocesshandler.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/clustercreator.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/lbgsampler.h \
 ../src/liftedinfer/liftedimportancesampling/includes/proposalconstructor.h \
 ../src/liftedinfer/liftedimportancesampling/includes/proposalstructure.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/filedump.h \
 ../src/liftedinfer/common/includes/stringconversionutils.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/ptpsearch.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/extensions.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/cache.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/unifier.h \
 ../src/liftedinfer/common/includes/lvrptpsearchtree.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/rulesutil.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/liftedimportancesampling/includes/lisapproxinference.h \
 ../src/liftedinfer/common/includes/propositionalresolution.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/extensions.h

../bin/obj/infer.o: ../src/infer/infer.cpp ../src/parser/fol.h ../src/logic/mln.h \
 ../src/logic/clause.h ../src/logic/domain.h ../src/util/equalstr.h \
 ../src/util/dualmap.h ../src/util/array.h ../src/util/random.h \
 ../src/util/equalstr.h ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicatetemplate.h ../src/util/array.h \
 ../src/logic/functiontemplate.h ../src/logic/predicate.h \
 ../src/logic/term.h ../src/util/hash.h ../src/util/hasharray.h \
 ../src/logic/function.h ../src/util/hashstring.h ../src/util/hasharray.h \
 ../src/util/hashlist.h ../src/util/arraysaccessor.h \
 ../src/util/powerset.h ../src/util/multdarray.h ../src/util/hashint.h \
 ../src/util/hash.h ../src/logic/database.h ../src/logic/groundpreds.h \
 ../src/infer/groundpredicate.h ../src/logic/domain.h \
 ../src/logic/predicate.h ../src/infer/groundclause.h \
 ../src/logic/clausesampler.h ../src/logic/truefalsegroundingsstore.h \
 ../src/util/random.h ../src/logic/clausehelper.h \
 ../src/logic/mlnhelper.h ../src/logic/database.h ../src/util/arguments.h \
 ../src/util/util.h ../src/infer/infer.h ../src/infer/mrf.h \
 ../src/util/timer.h ../src/infer/groundpredicate.h \
 ../src/learnwts/learnwts.h ../src/learnwts/indextranslator.h \
 ../src/infer/inferenceargs.h ../src/infer/inferutil.h \
 ../src/logic/clause.h ../src/infer/maxwalksat.h ../src/infer/sat.h \
 ../src/infer/inference.h ../src/logic/variablestate.h ../src/infer/mrf.h \
 ../src/logic/hvariablestate.h ../src/logic/Polynomial.h \
 ../src/learnwts/lbfgsp.h ../src/logic/Polynomial.h \
 ../src/logic/logsigmoid.h ../src/logic/logichelper.h \
 ../src/infer/maxwalksatparams.h ../src/infer/samplesatparams.h \
 ../src/infer/mcsat.h ../src/infer/mcmc.h ../src/infer/mcmcparams.h \
 ../src/infer/mcsatparams.h ../src/infer/unitpropagation.h \
 ../src/infer/gibbssampler.h ../src/infer/gibbsparams.h \
 ../src/infer/convergencetest.h ../src/util/meanvariance.h \
 ../src/infer/gelmanconvergencetest.h ../src/infer/simulatedtempering.h \
 ../src/infer/simulatedtemperingparams.h ../src/infer/bp.h \
 ../src/infer/bpparams.h ../src/infer/twowaymessage.h \
 ../src/infer/superclause.h ../src/infer/variable.h \
 ../src/infer/hypercube.h ../src/infer/hypercubeoperations.h \
 ../src/infer/hypercubereverseindex.h ../src/infer/hypercuberefinement.h \
 ../src/infer/auxfactor.h ../src/infer/factor.h ../src/infer/link.h \
 ../src/infer/node.h ../src/infer/superpred.h \
 ../src/infer/clausecounter.h ../src/infer/factorgraph.h \
 ../src/logic/clausefactory.h ../src/infer/hmcsat.h \
 ../src/infer/hmaxwalksat.h ../src/logic/logichelper.h

../bin/obj/learnstruct.o: ../src/learnstruct/learnstruct.cpp ../src/util/arguments.h \
 ../src/parser/fol.h ../src/logic/mln.h ../src/logic/clause.h \
 ../src/logic/domain.h ../src/util/equalstr.h ../src/util/dualmap.h \
 ../src/util/array.h ../src/util/random.h ../src/util/equalstr.h \
 ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicatetemplate.h ../src/util/array.h \
 ../src/logic/functiontemplate.h ../src/logic/predicate.h \
 ../src/logic/term.h ../src/util/hash.h ../src/util/hasharray.h \
 ../src/logic/function.h ../src/util/hashstring.h ../src/util/hasharray.h \
 ../src/util/hashlist.h ../src/util/arraysaccessor.h \
 ../src/util/powerset.h ../src/util/multdarray.h ../src/util/hashint.h \
 ../src/util/hash.h ../src/logic/database.h ../src/logic/groundpreds.h \
 ../src/infer/groundpredicate.h ../src/logic/domain.h \
 ../src/logic/predicate.h ../src/infer/groundclause.h \
 ../src/logic/clausesampler.h ../src/logic/truefalsegroundingsstore.h \
 ../src/util/random.h ../src/logic/clausehelper.h \
 ../src/logic/mlnhelper.h ../src/logic/database.h \
 ../src/learnwts/learnwts.h ../src/util/util.h ../src/util/timer.h \
 ../src/learnwts/indextranslator.h ../src/learnstruct/structlearn.h \
 ../src/logic/clausefactory.h ../src/learnwts/lbfgsb.h \
 ../src/learnwts/pseudologlikelihood.h ../src/logic/clause.h \
 ../src/logic/Polynomial.h ../src/learnwts/discriminativelearner.h \
 ../src/infer/infer.h ../src/infer/mrf.h ../src/infer/groundpredicate.h \
 ../src/infer/inferenceargs.h ../src/infer/inferutil.h \
 ../src/infer/maxwalksat.h ../src/infer/sat.h ../src/infer/inference.h \
 ../src/logic/variablestate.h ../src/infer/mrf.h \
 ../src/logic/hvariablestate.h ../src/logic/Polynomial.h \
 ../src/learnwts/lbfgsp.h ../src/logic/logsigmoid.h \
 ../src/logic/logichelper.h ../src/infer/maxwalksatparams.h \
 ../src/infer/samplesatparams.h ../src/infer/mcsat.h ../src/infer/mcmc.h \
 ../src/infer/mcmcparams.h ../src/infer/mcsatparams.h \
 ../src/infer/unitpropagation.h ../src/infer/gibbssampler.h \
 ../src/infer/gibbsparams.h ../src/infer/convergencetest.h \
 ../src/util/meanvariance.h ../src/infer/gelmanconvergencetest.h \
 ../src/infer/simulatedtempering.h \
 ../src/infer/simulatedtemperingparams.h ../src/infer/bp.h \
 ../src/infer/bpparams.h ../src/infer/twowaymessage.h \
 ../src/infer/superclause.h ../src/infer/variable.h \
 ../src/infer/hypercube.h ../src/infer/hypercubeoperations.h \
 ../src/infer/hypercubereverseindex.h ../src/infer/hypercuberefinement.h \
 ../src/infer/auxfactor.h ../src/infer/factor.h ../src/infer/link.h \
 ../src/infer/node.h ../src/infer/superpred.h \
 ../src/infer/clausecounter.h ../src/infer/factorgraph.h \
 ../src/infer/hmcsat.h ../src/infer/hmaxwalksat.h \
 ../src/logic/logichelper.h ../src/infer/maxwalksat.h ../src/infer/bp.h

../bin/obj/Polynomial.o: ../src/logic/Polynomial.cpp ../src/logic/Polynomial.h \
 ../src/util/array.h ../src/util/random.h ../src/util/complex.h \
 ../src/util/array.h

../bin/obj/logichelper.o: ../src/logic/logichelper.cpp ../src/logic/logichelper.h \
 ../src/util/array.h ../src/util/random.h

../bin/obj/inferutil.o: ../src/infer/inferutil.cpp ../src/infer/inferutil.h \
 ../src/util/array.h ../src/util/random.h ../src/util/hashint.h \
 ../src/util/hasharray.h ../src/util/hash.h ../src/util/array.h \
 ../src/logic/clause.h ../src/logic/domain.h ../src/util/equalstr.h \
 ../src/util/dualmap.h ../src/util/equalstr.h ../src/util/constdualmap.h \
 ../src/util/strint.h ../src/logic/predicatetemplate.h \
 ../src/logic/functiontemplate.h ../src/logic/predicate.h \
 ../src/logic/term.h ../src/util/hash.h ../src/util/hasharray.h \
 ../src/logic/function.h ../src/util/hashstring.h ../src/util/hashlist.h \
 ../src/util/arraysaccessor.h ../src/util/powerset.h \
 ../src/util/multdarray.h ../src/logic/database.h \
 ../src/logic/groundpreds.h ../src/infer/groundpredicate.h \
 ../src/logic/domain.h ../src/logic/predicate.h \
 ../src/infer/groundclause.h ../src/logic/clausesampler.h \
 ../src/logic/truefalsegroundingsstore.h ../src/util/random.h \
 ../src/logic/clausehelper.h

../bin/obj/createhypercube.o: ../src/infer/createhypercube.cpp \
 ../src/infer/hypercube.h ../src/util/util.h ../src/util/array.h \
 ../src/util/random.h ../src/util/hashint.h ../src/util/hasharray.h \
 ../src/util/hash.h ../src/util/array.h ../src/infer/inferutil.h \
 ../src/logic/clause.h ../src/logic/domain.h ../src/util/equalstr.h \
 ../src/util/dualmap.h ../src/util/equalstr.h ../src/util/constdualmap.h \
 ../src/util/strint.h ../src/logic/predicatetemplate.h \
 ../src/logic/functiontemplate.h ../src/logic/predicate.h \
 ../src/logic/term.h ../src/util/hash.h ../src/util/hasharray.h \
 ../src/logic/function.h ../src/util/hashstring.h ../src/util/hashlist.h \
 ../src/util/arraysaccessor.h ../src/util/powerset.h \
 ../src/util/multdarray.h ../src/logic/database.h \
 ../src/logic/groundpreds.h ../src/infer/groundpredicate.h \
 ../src/logic/domain.h ../src/logic/predicate.h \
 ../src/infer/groundclause.h ../src/logic/clausesampler.h \
 ../src/logic/truefalsegroundingsstore.h ../src/util/random.h \
 ../src/logic/clausehelper.h ../src/infer/hypercuberefinement.h \
 ../src/logic/mln.h ../src/logic/clause.h ../src/logic/mlnhelper.h \
 ../src/infer/hypercubeoperations.h ../src/infer/hypercubereverseindex.h \
 ../src/infer/mrf.h ../src/util/timer.h ../src/infer/groundpredicate.h

../bin/obj/createhypercubebasic.o: ../src/infer/createhypercubebasic.cpp \
 ../src/infer/hypercube.h ../src/util/util.h ../src/util/array.h \
 ../src/util/random.h ../src/util/hashint.h ../src/util/hasharray.h \
 ../src/util/hash.h ../src/util/array.h ../src/infer/inferutil.h \
 ../src/logic/clause.h ../src/logic/domain.h ../src/util/equalstr.h \
 ../src/util/dualmap.h ../src/util/equalstr.h ../src/util/constdualmap.h \
 ../src/util/strint.h ../src/logic/predicatetemplate.h \
 ../src/logic/functiontemplate.h ../src/logic/predicate.h \
 ../src/logic/term.h ../src/util/hash.h ../src/util/hasharray.h \
 ../src/logic/function.h ../src/util/hashstring.h ../src/util/hashlist.h \
 ../src/util/arraysaccessor.h ../src/util/powerset.h \
 ../src/util/multdarray.h ../src/logic/database.h \
 ../src/logic/groundpreds.h ../src/infer/groundpredicate.h \
 ../src/logic/domain.h ../src/logic/predicate.h \
 ../src/infer/groundclause.h ../src/logic/clausesampler.h \
 ../src/logic/truefalsegroundingsstore.h ../src/util/random.h \
 ../src/logic/clausehelper.h ../src/infer/hypercuberefinement.h \
 ../src/logic/mln.h ../src/logic/clause.h ../src/logic/mlnhelper.h

../bin/obj/createhypercubedt.o: ../src/infer/createhypercubedt.cpp \
 ../src/util/hashint.h ../src/util/hasharray.h ../src/util/random.h \
 ../src/util/hash.h ../src/util/array.h ../src/infer/inferutil.h \
 ../src/util/array.h ../src/logic/clause.h ../src/logic/domain.h \
 ../src/util/equalstr.h ../src/util/dualmap.h ../src/util/equalstr.h \
 ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicatetemplate.h ../src/logic/functiontemplate.h \
 ../src/logic/predicate.h ../src/logic/term.h ../src/util/hash.h \
 ../src/util/hasharray.h ../src/logic/function.h ../src/util/hashstring.h \
 ../src/util/hashlist.h ../src/util/arraysaccessor.h \
 ../src/util/powerset.h ../src/util/multdarray.h ../src/logic/database.h \
 ../src/logic/groundpreds.h ../src/infer/groundpredicate.h \
 ../src/logic/domain.h ../src/logic/predicate.h \
 ../src/infer/groundclause.h ../src/logic/clausesampler.h \
 ../src/logic/truefalsegroundingsstore.h ../src/util/random.h \
 ../src/logic/clausehelper.h ../src/infer/hypercube.h ../src/util/util.h \
 ../src/infer/hypercuberefinement.h ../src/logic/mln.h \
 ../src/logic/clause.h ../src/logic/mlnhelper.h

../bin/obj/joinhypercube.o: ../src/infer/joinhypercube.cpp ../src/util/util.h \
 ../src/infer/hypercube.h ../src/util/array.h ../src/util/random.h \
 ../src/util/hashint.h ../src/util/hasharray.h ../src/util/hash.h \
 ../src/util/array.h ../src/infer/inferutil.h ../src/logic/clause.h \
 ../src/logic/domain.h ../src/util/equalstr.h ../src/util/dualmap.h \
 ../src/util/equalstr.h ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicatetemplate.h ../src/logic/functiontemplate.h \
 ../src/logic/predicate.h ../src/logic/term.h ../src/util/hash.h \
 ../src/util/hasharray.h ../src/logic/function.h ../src/util/hashstring.h \
 ../src/util/hashlist.h ../src/util/arraysaccessor.h \
 ../src/util/powerset.h ../src/util/multdarray.h ../src/logic/database.h \
 ../src/logic/groundpreds.h ../src/infer/groundpredicate.h \
 ../src/logic/domain.h ../src/logic/predicate.h \
 ../src/infer/groundclause.h ../src/logic/clausesampler.h \
 ../src/logic/truefalsegroundingsstore.h ../src/util/random.h \
 ../src/logic/clausehelper.h ../src/infer/hypercubereverseindex.h \
 ../src/infer/mrf.h ../src/util/timer.h ../src/logic/mln.h \
 ../src/logic/clause.h ../src/logic/mlnhelper.h \
 ../src/infer/groundpredicate.h ../src/infer/partialgndhctuple.h \
 ../src/infer/hypercubeoperations.h ../src/infer/hypercuberefinement.h \
 ../src/infer/superclause.h ../src/infer/variable.h

../bin/obj/refinehypercube.o: ../src/infer/refinehypercube.cpp ../src/util/util.h \
 ../src/infer/hypercube.h ../src/util/array.h ../src/util/random.h \
 ../src/util/hashint.h ../src/util/hasharray.h ../src/util/hash.h \
 ../src/util/array.h ../src/infer/inferutil.h ../src/logic/clause.h \
 ../src/logic/domain.h ../src/util/equalstr.h ../src/util/dualmap.h \
 ../src/util/equalstr.h ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicatetemplate.h ../src/logic/functiontemplate.h \
 ../src/logic/predicate.h ../src/logic/term.h ../src/util/hash.h \
 ../src/util/hasharray.h ../src/logic/function.h ../src/util/hashstring.h \
 ../src/util/hashlist.h ../src/util/arraysaccessor.h \
 ../src/util/powerset.h ../src/util/multdarray.h ../src/logic/database.h \
 ../src/logic/groundpreds.h ../src/infer/groundpredicate.h \
 ../src/logic/domain.h ../src/logic/predicate.h \
 ../src/infer/groundclause.h ../src/logic/clausesampler.h \
 ../src/logic/truefalsegroundingsstore.h ../src/util/random.h \
 ../src/logic/clausehelper.h ../src/infer/hypercuberefinement.h \
 ../src/logic/mln.h ../src/logic/clause.h ../src/logic/mlnhelper.h

../bin/obj/superpred.o: ../src/infer/superpred.cpp ../src/util/util.h \
 ../src/infer/mrf.h ../src/util/timer.h ../src/logic/mln.h \
 ../src/logic/clause.h ../src/logic/domain.h ../src/util/equalstr.h \
 ../src/util/dualmap.h ../src/util/array.h ../src/util/random.h \
 ../src/util/equalstr.h ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicatetemplate.h ../src/util/array.h \
 ../src/logic/functiontemplate.h ../src/logic/predicate.h \
 ../src/logic/term.h ../src/util/hash.h ../src/util/hasharray.h \
 ../src/logic/function.h ../src/util/hashstring.h ../src/util/hasharray.h \
 ../src/util/hashlist.h ../src/util/arraysaccessor.h \
 ../src/util/powerset.h ../src/util/multdarray.h ../src/util/hashint.h \
 ../src/util/hash.h ../src/logic/database.h ../src/logic/groundpreds.h \
 ../src/infer/groundpredicate.h ../src/logic/domain.h \
 ../src/logic/predicate.h ../src/infer/groundclause.h \
 ../src/logic/clausesampler.h ../src/logic/truefalsegroundingsstore.h \
 ../src/util/random.h ../src/logic/clausehelper.h \
 ../src/logic/mlnhelper.h ../src/infer/groundpredicate.h \
 ../src/infer/superclause.h ../src/infer/variable.h \
 ../src/infer/hypercube.h ../src/infer/inferutil.h ../src/logic/clause.h \
 ../src/infer/hypercubeoperations.h ../src/infer/hypercubereverseindex.h \
 ../src/infer/hypercuberefinement.h ../src/infer/superpred.h \
 ../src/infer/clausecounter.h

../bin/obj/superclause.o: ../src/infer/superclause.cpp ../src/util/util.h \
 ../src/infer/mrf.h ../src/util/timer.h ../src/logic/mln.h \
 ../src/logic/clause.h ../src/logic/domain.h ../src/util/equalstr.h \
 ../src/util/dualmap.h ../src/util/array.h ../src/util/random.h \
 ../src/util/equalstr.h ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicatetemplate.h ../src/util/array.h \
 ../src/logic/functiontemplate.h ../src/logic/predicate.h \
 ../src/logic/term.h ../src/util/hash.h ../src/util/hasharray.h \
 ../src/logic/function.h ../src/util/hashstring.h ../src/util/hasharray.h \
 ../src/util/hashlist.h ../src/util/arraysaccessor.h \
 ../src/util/powerset.h ../src/util/multdarray.h ../src/util/hashint.h \
 ../src/util/hash.h ../src/logic/database.h ../src/logic/groundpreds.h \
 ../src/infer/groundpredicate.h ../src/logic/domain.h \
 ../src/logic/predicate.h ../src/infer/groundclause.h \
 ../src/logic/clausesampler.h ../src/logic/truefalsegroundingsstore.h \
 ../src/util/random.h ../src/logic/clausehelper.h \
 ../src/logic/mlnhelper.h ../src/infer/groundpredicate.h \
 ../src/infer/superclause.h ../src/infer/variable.h \
 ../src/infer/hypercube.h ../src/infer/inferutil.h ../src/logic/clause.h \
 ../src/infer/hypercubeoperations.h ../src/infer/hypercubereverseindex.h \
 ../src/infer/hypercuberefinement.h ../src/infer/superpred.h \
 ../src/infer/clausecounter.h

../bin/obj/node.o: ../src/infer/node.cpp ../src/infer/node.h ../src/util/util.h \
 ../src/infer/mrf.h ../src/util/timer.h ../src/logic/mln.h \
 ../src/logic/clause.h ../src/logic/domain.h ../src/util/equalstr.h \
 ../src/util/dualmap.h ../src/util/array.h ../src/util/random.h \
 ../src/util/equalstr.h ../src/util/constdualmap.h ../src/util/strint.h \
 ../src/logic/predicatetemplate.h ../src/util/array.h \
 ../src/logic/functiontemplate.h ../src/logic/predicate.h \
 ../src/logic/term.h ../src/util/hash.h ../src/util/hasharray.h \
 ../src/logic/function.h ../src/util/hashstring.h ../src/util/hasharray.h \
 ../src/util/hashlist.h ../src/util/arraysaccessor.h \
 ../src/util/powerset.h ../src/util/multdarray.h ../src/util/hashint.h \
 ../src/util/hash.h ../src/logic/database.h ../src/logic/groundpreds.h \
 ../src/infer/groundpredicate.h ../src/logic/domain.h \
 ../src/logic/predicate.h ../src/infer/groundclause.h \
 ../src/logic/clausesampler.h ../src/logic/truefalsegroundingsstore.h \
 ../src/util/random.h ../src/logic/clausehelper.h \
 ../src/logic/mlnhelper.h ../src/infer/groundpredicate.h \
 ../src/infer/link.h ../src/infer/superpred.h \
 ../src/infer/clausecounter.h ../src/infer/inferutil.h \
 ../src/logic/clause.h ../src/infer/superclause.h ../src/infer/variable.h \
 ../src/infer/hypercube.h ../src/infer/hypercubeoperations.h \
 ../src/infer/hypercubereverseindex.h ../src/infer/hypercuberefinement.h \
 ../src/infer/twowaymessage.h ../src/infer/factor.h

../bin/obj/factor.o: ../src/infer/factor.cpp ../src/infer/factor.h \
 ../src/util/util.h ../src/infer/mrf.h ../src/util/timer.h \
 ../src/logic/mln.h ../src/logic/clause.h ../src/logic/domain.h \
 ../src/util/equalstr.h ../src/util/dualmap.h ../src/util/array.h \
 ../src/util/random.h ../src/util/equalstr.h ../src/util/constdualmap.h \
 ../src/util/strint.h ../src/logic/predicatetemplate.h \
 ../src/util/array.h ../src/logic/functiontemplate.h \
 ../src/logic/predicate.h ../src/logic/term.h ../src/util/hash.h \
 ../src/util/hasharray.h ../src/logic/function.h ../src/util/hashstring.h \
 ../src/util/hasharray.h ../src/util/hashlist.h \
 ../src/util/arraysaccessor.h ../src/util/powerset.h \
 ../src/util/multdarray.h ../src/util/hashint.h ../src/util/hash.h \
 ../src/logic/database.h ../src/logic/groundpreds.h \
 ../src/infer/groundpredicate.h ../src/logic/domain.h \
 ../src/logic/predicate.h ../src/infer/groundclause.h \
 ../src/logic/clausesampler.h ../src/logic/truefalsegroundingsstore.h \
 ../src/util/random.h ../src/logic/clausehelper.h \
 ../src/logic/mlnhelper.h ../src/infer/groundpredicate.h \
 ../src/infer/link.h ../src/infer/node.h ../src/infer/superpred.h \
 ../src/infer/clausecounter.h ../src/infer/inferutil.h \
 ../src/logic/clause.h ../src/infer/superclause.h ../src/infer/variable.h \
 ../src/infer/hypercube.h ../src/infer/hypercubeoperations.h \
 ../src/infer/hypercubereverseindex.h ../src/infer/hypercuberefinement.h \
 ../src/infer/twowaymessage.h

../bin/obj/liftedalgsconvertor.o: \
 ../src/liftedinfer/bridge/src/liftedalgsconvertor.cpp \
 ../src/liftedinfer/bridge/includes/liftedalgsconvertor.h \
 ../src/logic/mln.h ../src/logic/clause.h ../src/logic/domain.h \
 ../src/util/equalstr.h ../src/util/dualmap.h ../src/util/array.h \
 ../src/util/random.h ../src/util/equalstr.h ../src/util/constdualmap.h \
 ../src/util/strint.h ../src/logic/predicatetemplate.h \
 ../src/util/array.h ../src/logic/functiontemplate.h \
 ../src/logic/predicate.h ../src/logic/term.h ../src/util/hash.h \
 ../src/util/hasharray.h ../src/logic/function.h ../src/util/hashstring.h \
 ../src/util/hasharray.h ../src/util/hashlist.h \
 ../src/util/arraysaccessor.h ../src/util/powerset.h \
 ../src/util/multdarray.h ../src/util/hashint.h ../src/util/hash.h \
 ../src/logic/database.h ../src/logic/groundpreds.h \
 ../src/infer/groundpredicate.h ../src/logic/domain.h \
 ../src/logic/predicate.h ../src/infer/groundclause.h \
 ../src/logic/clausesampler.h ../src/logic/truefalsegroundingsstore.h \
 ../src/util/random.h ../src/logic/clausehelper.h \
 ../src/logic/mlnhelper.h ../src/logic/variablestate.h ../src/infer/mrf.h \
 ../src/util/timer.h ../src/infer/groundpredicate.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/bridge/includes/liftedalgshandler.h \
 ../src/liftedinfer/common/includes/lvparams.h ../src/infer/mcmcparams.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/clustercreator.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/ptpsampler.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/blockexactinference.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/common/includes/heuristics.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/lvrnormpropagation.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/randomgenutil.h \
 ../src/liftedinfer/common/includes/lvrsingletonnormpropagation.h \
 ../src/liftedinfer/common/includes/lvrptptreesampling.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/lvrcluster.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/gibbsprocesshandler.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/clustercreator.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/lbgsampler.h \
 ../src/liftedinfer/liftedimportancesampling/includes/proposalconstructor.h \
 ../src/liftedinfer/liftedimportancesampling/includes/proposalstructure.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/filedump.h \
 ../src/liftedinfer/common/includes/stringconversionutils.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/ptpsearch.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/extensions.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/cache.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/unifier.h \
 ../src/liftedinfer/common/includes/lvrptpsearchtree.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/rulesutil.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/liftedimportancesampling/includes/lisapproxinference.h \
 ../src/liftedinfer/common/includes/propositionalresolution.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/extensions.h \
 ../src/liftedinfer/common/includes/queryupdater.h

../bin/obj/liftedalgshandler.o: ../src/liftedinfer/bridge/src/liftedalgshandler.cpp \
 ../src/liftedinfer/bridge/includes/liftedalgshandler.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/lvparams.h ../src/infer/mcmcparams.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/clustercreator.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/ptpsampler.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/blockexactinference.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/common/includes/heuristics.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/lvrnormpropagation.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/randomgenutil.h \
 ../src/liftedinfer/common/includes/lvrsingletonnormpropagation.h \
 ../src/liftedinfer/common/includes/lvrptptreesampling.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/lvrcluster.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/gibbsprocesshandler.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/clustercreator.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/lbgsampler.h \
 ../src/liftedinfer/liftedimportancesampling/includes/proposalconstructor.h \
 ../src/liftedinfer/liftedimportancesampling/includes/proposalstructure.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/filedump.h \
 ../src/liftedinfer/common/includes/stringconversionutils.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/ptpsearch.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/extensions.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/cache.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/unifier.h \
 ../src/liftedinfer/common/includes/lvrptpsearchtree.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/rulesutil.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/liftedimportancesampling/includes/lisapproxinference.h \
 ../src/liftedinfer/common/includes/propositionalresolution.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/extensions.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/wmconvertor.h

../bin/obj/decomposer.o: ../src/liftedinfer/common/src/decomposer.cpp \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/cleanuputils.h

../bin/obj/fileutils.o: ../src/liftedinfer/common/src/fileutils.cpp \
 ../src/liftedinfer/common/includes/fileutils.h \
 ../src/liftedinfer/common/includes/logdouble.h

../bin/obj/heuristics.o: ../src/liftedinfer/common/src/heuristics.cpp \
 ../src/liftedinfer/common/includes/heuristics.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h

../bin/obj/logdouble.o: ../src/liftedinfer/common/src/logdouble.cpp \
 ../src/liftedinfer/common/includes/logdouble.h

../bin/obj/lvrmln.o: ../src/liftedinfer/common/src/lvrmln.cpp \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h

../bin/obj/mathutils.o: ../src/liftedinfer/common/src/mathutils.cpp \
 ../src/liftedinfer/common/includes/mathutils.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/randomgenutil.h

../bin/obj/normalizer.o: ../src/liftedinfer/common/src/normalizer.cpp \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/src/setutils.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/rulesutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/normpreservinggrounder.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/randomgenutil.h

../bin/obj/propositionalresolution.o: \
 ../src/liftedinfer/common/src/propositionalresolution.cpp \
 ../src/liftedinfer/common/includes/propositionalresolution.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/cleanuputils.h

../bin/obj/queryupdater.o: ../src/liftedinfer/common/src/queryupdater.cpp \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/randomgenutil.h \
 ../src/liftedinfer/common/includes/fileutils.h

../bin/obj/samplingalgs.o: ../src/liftedinfer/common/src/samplingalgs.cpp \
 ../src/liftedinfer/common/includes/samplingalgs.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedimportancesampling/includes/proposalstructure.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/lvparams.h ../src/infer/mcmcparams.h \
 ../src/liftedinfer/common/includes/lvparams.h \
 ../src/liftedinfer/common/includes/mathutils.h \
 ../src/liftedinfer/common/includes/randomgenutil.h

../bin/obj/randomgenutil.o: ../src/liftedinfer/common/src/randomgenutil.cpp \
 ../src/liftedinfer/common/includes/randomgenutil.h

../bin/obj/lvrptptreesampling.o: \
 ../src/liftedinfer/common/src/lvrptptreesampling.cpp \
 ../src/liftedinfer/common/includes/lvrptptreesampling.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/randomgenutil.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h

../bin/obj/lvrptpsearchtree.o: ../src/liftedinfer/common/src/lvrptpsearchtree.cpp \
 ../src/liftedinfer/common/includes/lvrptpsearchtree.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/mathutils.h

../bin/obj/blockexactinference.o: \
 ../src/liftedinfer/liftedblockedgibbs/src/blockexactinference.cpp \
 ../src/liftedinfer/liftedblockedgibbs/includes/blockexactinference.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/common/includes/heuristics.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/lvrnormpropagation.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/randomgenutil.h \
 ../src/liftedinfer/common/includes/lvrsingletonnormpropagation.h \
 ../src/liftedinfer/common/includes/lvrptptreesampling.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/propositionalresolution.h \
 ../src/liftedinfer/common/includes/randomgenutil.h

../bin/obj/clustercreator.o: \
 ../src/liftedinfer/liftedblockedgibbs/src/clustercreator.cpp \
 ../src/liftedinfer/liftedblockedgibbs/includes/clustercreator.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/ptpsampler.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/blockexactinference.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/common/includes/heuristics.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/lvrnormpropagation.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/randomgenutil.h \
 ../src/liftedinfer/common/includes/lvrsingletonnormpropagation.h \
 ../src/liftedinfer/common/includes/lvrptptreesampling.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/lvrcluster.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/infer/mcmcparams.h ../src/liftedinfer/common/includes/lvparams.h \
 ../src/liftedinfer/common/includes/randomgenutil.h

../bin/obj/lbgsampler.o: ../src/liftedinfer/liftedblockedgibbs/src/lbgsampler.cpp \
 ../src/liftedinfer/liftedblockedgibbs/includes/lbgsampler.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/ptpsampler.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/blockexactinference.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/common/includes/heuristics.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/lvrnormpropagation.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/randomgenutil.h \
 ../src/liftedinfer/common/includes/lvrsingletonnormpropagation.h \
 ../src/liftedinfer/common/includes/lvrptptreesampling.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/lvrcluster.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/lvparams.h ../src/infer/mcmcparams.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/stringconversionutils.h

../bin/obj/lvrcluster.o: ../src/liftedinfer/liftedblockedgibbs/src/lvrcluster.cpp \
 ../src/liftedinfer/liftedblockedgibbs/includes/lvrcluster.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/randomgenutil.h

../bin/obj/ptpsampler.o: ../src/liftedinfer/liftedblockedgibbs/src/ptpsampler.cpp \
 ../src/liftedinfer/liftedblockedgibbs/includes/ptpsampler.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/blockexactinference.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/common/includes/heuristics.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/lvrnormpropagation.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/randomgenutil.h \
 ../src/liftedinfer/common/includes/lvrsingletonnormpropagation.h \
 ../src/liftedinfer/common/includes/lvrptptreesampling.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/lvrcluster.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h

../bin/obj/lisapproxinference.o: \
 ../src/liftedinfer/liftedimportancesampling/src/lisapproxinference.cpp \
 ../src/liftedinfer/liftedimportancesampling/includes/lisapproxinference.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/propositionalresolution.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/extensions.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/cache.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/unifier.h \
 ../src/liftedinfer/common/includes/heuristics.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/liftedimportancesampling/includes/proposalstructure.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvparams.h ../src/infer/mcmcparams.h \
 ../src/liftedinfer/common/includes/lvrnormpropagation.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/randomgenutil.h \
 ../src/liftedinfer/common/includes/lvrsingletonnormpropagation.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/rulesutil.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/mathutils.h \
 ../src/liftedinfer/common/includes/samplingalgs.h \
 ../src/liftedinfer/liftedimportancesampling/includes/proposalstructure.h \
 ../src/liftedinfer/common/includes/lvparams.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h

../bin/obj/proposalconstructor.o: \
 ../src/liftedinfer/liftedimportancesampling/src/proposalconstructor.cpp \
 ../src/liftedinfer/liftedimportancesampling/includes/proposalconstructor.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/liftedimportancesampling/includes/proposalstructure.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/lvparams.h ../src/infer/mcmcparams.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/filedump.h \
 ../src/liftedinfer/common/includes/stringconversionutils.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/ptpsearch.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/extensions.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/cache.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/unifier.h \
 ../src/liftedinfer/common/includes/heuristics.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrptptreesampling.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/lvrptpsearchtree.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/rulesutil.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/liftedimportancesampling/includes/lisapproxinference.h \
 ../src/liftedinfer/common/includes/propositionalresolution.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/extensions.h \
 ../src/liftedinfer/common/includes/lvrnormpropagation.h \
 ../src/liftedinfer/common/includes/randomgenutil.h \
 ../src/liftedinfer/common/includes/lvrsingletonnormpropagation.h

../bin/obj/proposalstructure.o: \
 ../src/liftedinfer/liftedimportancesampling/src/proposalstructure.cpp \
 ../src/liftedinfer/liftedimportancesampling/includes/proposalstructure.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/lvparams.h ../src/infer/mcmcparams.h \
 ../src/liftedinfer/common/includes/filedump.h \
 ../src/liftedinfer/common/includes/stringconversionutils.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/randomgenutil.h

../bin/obj/cache.o: ../src/liftedinfer/liftedweightedmodelcounting/src/cache.cpp \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/cache.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/unifier.h \
 ../src/liftedinfer/common/includes/cleanuputils.h

../bin/obj/extensions.o: \
 ../src/liftedinfer/liftedweightedmodelcounting/src/extensions.cpp \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/extensions.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/cache.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/unifier.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h

../bin/obj/ptpsearch.o: \
 ../src/liftedinfer/liftedweightedmodelcounting/src/ptpsearch.cpp \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/ptpsearch.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/extensions.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/cache.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/unifier.h \
 ../src/liftedinfer/common/includes/heuristics.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/lvparams.h ../src/infer/mcmcparams.h \
 ../src/liftedinfer/common/includes/lvrptptreesampling.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/lvrptpsearchtree.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/propositionalresolution.h \
 ../src/liftedinfer/common/includes/rulesutil.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvrsingletonnormpropagation.h

../bin/obj/resolver.o: \
 ../src/liftedinfer/liftedweightedmodelcounting/src/resolver.cpp \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h

../bin/obj/sampler.o: ../src/liftedinfer/liftedweightedmodelcounting/src/sampler.cpp \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/randomgenutil.h

../bin/obj/splitter.o: \
 ../src/liftedinfer/liftedweightedmodelcounting/src/splitter.cpp \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/common/includes/cleanuputils.h

../bin/obj/unifier.o: ../src/liftedinfer/liftedweightedmodelcounting/src/unifier.cpp \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/unifier.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/cleanuputils.h

../bin/obj/wmconvertor.o: \
 ../src/liftedinfer/liftedweightedmodelcounting/src/wmconvertor.cpp \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/wmconvertor.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h

../bin/obj/parser.o: ../src/liftedinfer/tests/parser.cpp \
 ../src/liftedinfer/tests/parser.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/stringconversionutils.h

../bin/obj/testmodelcount.o: ../src/liftedinfer/tests/testmodelcount.cpp \
 ../src/liftedinfer/tests/testmodelcount.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/ptpsearch.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/extensions.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/cache.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/unifier.h \
 ../src/liftedinfer/common/includes/heuristics.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/lvparams.h ../src/infer/mcmcparams.h \
 ../src/liftedinfer/common/includes/lvrptptreesampling.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/lvrptpsearchtree.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/wmconvertor.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/tests/parser.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/blockexactinference.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/lvrnormpropagation.h \
 ../src/liftedinfer/common/includes/randomgenutil.h \
 ../src/liftedinfer/common/includes/lvrsingletonnormpropagation.h

../bin/obj/testformulainfer.o: ../src/liftedinfer/tests/testformulainfer.cpp \
 ../src/liftedinfer/tests/testformulainfer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/common/includes/paramconstants.h \
 ../src/liftedinfer/liftedblockedgibbs/includes/blockexactinference.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/lvrmln.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/sampler.h \
 ../src/liftedinfer/common/includes/logdouble.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/resolver.h \
 ../src/liftedinfer/common/includes/heuristics.h \
 ../src/liftedinfer/common/includes/decomposer.h \
 ../src/liftedinfer/common/includes/normalizer.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/clusterutil.h \
 ../src/liftedinfer/common/includes/cleanuputils.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/common/includes/lvrnormpropagation.h \
 ../src/liftedinfer/common/includes/hashalgorithm.h \
 ../src/liftedinfer/common/includes/lvrpermutations.h \
 ../src/liftedinfer/common/includes/randomgenutil.h \
 ../src/liftedinfer/common/includes/lvrsingletonnormpropagation.h \
 ../src/liftedinfer/common/includes/lvrptptreesampling.h \
 ../src/liftedinfer/common/includes/ptptree.h \
 ../src/liftedinfer/tests/parser.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/ptpsearch.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/splitter.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/extensions.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/cache.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/unifier.h \
 ../src/liftedinfer/common/includes/lvparams.h ../src/infer/mcmcparams.h \
 ../src/liftedinfer/common/includes/lvrptpsearchtree.h \
 ../src/liftedinfer/common/includes/queryupdater.h \
 ../src/liftedinfer/common/includes/lvratomhashtemplate.h \
 ../src/liftedinfer/liftedweightedmodelcounting/includes/wmconvertor.h \
 ../src/liftedinfer/common/includes/stringconversionutils.h

