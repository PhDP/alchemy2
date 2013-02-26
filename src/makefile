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
