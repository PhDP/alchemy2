/*
 * All of the documentation and software included in the
 * Alchemy Software is copyrighted by Stanley Kok, Parag
 * Singla, Matthew Richardson, Pedro Domingos, Marc
 * Sumner, Hoifung Poon, Daniel Lowd, and Jue Wang.
 * 
 * Copyright [2004-09] Stanley Kok, Parag Singla, Matthew
 * Richardson, Pedro Domingos, Marc Sumner, Hoifung
 * Poon, Daniel Lowd, and Jue Wang. All rights reserved.
 * 
 * Contact: Pedro Domingos, University of Washington
 * (pedrod@cs.washington.edu).
 * 
 * Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that
 * the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above
 * copyright notice, this list of conditions and the
 * following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the
 * above copyright notice, this list of conditions and the
 * following disclaimer in the documentation and/or other
 * materials provided with the distribution.
 * 
 * 3. All advertising materials mentioning features or use
 * of this software must display the following
 * acknowledgment: "This product includes software
 * developed by Stanley Kok, Parag Singla, Matthew
 * Richardson, Pedro Domingos, Marc Sumner, Hoifung
 * Poon, Daniel Lowd, and Jue Wang in the Department of
 * Computer Science and Engineering at the University of
 * Washington".
 * 
 * 4. Your publications acknowledge the use or
 * contribution made by the Software to your research
 * using the following citation(s): 
 * Stanley Kok, Parag Singla, Matthew Richardson and
 * Pedro Domingos (2005). "The Alchemy System for
 * Statistical Relational AI", Technical Report,
 * Department of Computer Science and Engineering,
 * University of Washington, Seattle, WA.
 * http://alchemy.cs.washington.edu.
 * 
 * 5. Neither the name of the University of Washington nor
 * the names of its contributors may be used to endorse or
 * promote products derived from this software without
 * specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF WASHINGTON
 * AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OF WASHINGTON OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
#ifndef DOMAIN_H_JUN_21_2005
#define DOMAIN_H_JUN_21_2005

#include <ext/hash_set>
#include "equalstr.h"
#include "dualmap.h"
#include "constdualmap.h"
#include "predicatetemplate.h"
#include "functiontemplate.h"
#include "predicate.h"
#include "function.h"

class Clause;
class Database;
class TrueFalseGroundingsStore;
class MLN;


typedef hash_map<const char*, const PredicateTemplate*, hash<const char*>,
                 EqualStr>  StrToPredTemplateMap;

typedef hash_map<const char*, const FunctionTemplate*, hash<const char*>,
                 EqualStr>  StrToFuncTemplateMap;


class Domain
{
 public:
    // Name of type used in a unary predicate to simulate propositional constant
  static const char* PROPOSITIONAL_TYPE;
    // Name of constant used in a unary predicate to simulate propositional
    // constant
  static const char* PROPOSITIONAL_CONSTANT;

 public:
  Domain() : typeDualMap_(new DualMap), constDualMap_(new ConstDualMap),
             predDualMap_(new DualMap), funcDualMap_(new DualMap),
             strToPredTemplateMap_(new StrToPredTemplateMap),
             strToFuncTemplateMap_(new StrToFuncTemplateMap),
             constantsByType_(new Array<Array<int>*>),
             externalConstantsByType_(new Array<Array<int>*>), db_(NULL),
             trueFalseGroundingsStore_(NULL), funcSet_(new FunctionSet),
             predBlocks_(new Array<Predicate*>),
             truePredsInBlock_(new Array<Predicate*>),
             blockSizes_(new Array<int>),
             blockEvidence_(new Array<bool>),
             //externalConstant_(new Array<bool>),
             numNonEvidAtomsPerPred_(new Array<int>),
             numTrueNonEvidGndingsPerClause_(new Array<double>),             
             numFalseNonEvidGndingsPerClause_(new Array<double>)             
  {
    equalPredTemplate_ = new PredicateTemplate();
    equalPredTemplate_->setName(PredicateTemplate::EQUAL_NAME);
    emptyPredTemplate_ = new PredicateTemplate();
    emptyPredTemplate_->setName(PredicateTemplate::EMPTY_NAME);
    emptyFuncUnaryTemplate_ = new FunctionTemplate();
    emptyFuncUnaryTemplate_->setName(FunctionTemplate::EMPTY_FTEMPLATE_NAME);
    emptyFuncBinaryTemplate_ = new FunctionTemplate();
    emptyFuncBinaryTemplate_->setName(FunctionTemplate::EMPTY_FTEMPLATE_NAME);
    
    int typeId = addType(PredicateTemplate::ANY_TYPE_NAME);
    if (typeId != 0)
    {
      cout << "ERROR in domain.h: typeId of EQUAL predicate is not 0." << endl;
      exit(-1);
    }
    equalPredTemplate_->appendTermType(PredicateTemplate::ANY_TYPE_NAME, 
                                       false, this);
    equalPredTemplate_->appendTermType(PredicateTemplate::ANY_TYPE_NAME, 
                                       false, this);
    emptyPredTemplate_->appendTermType(PredicateTemplate::ANY_TYPE_NAME, 
                                       false, this);
    emptyPredTemplate_->appendTermType(PredicateTemplate::ANY_TYPE_NAME, 
                                       false, this);
    emptyFuncUnaryTemplate_->appendTermType(FunctionTemplate::ANY_TYPE_NAME, 
                                       false, this);
    emptyFuncBinaryTemplate_->appendTermType(FunctionTemplate::ANY_TYPE_NAME, 
                                       false, this);
    emptyFuncBinaryTemplate_->appendTermType(FunctionTemplate::ANY_TYPE_NAME, 
                                       false, this);                                       
  }


  ~Domain();


  void replaceTypeDualMap(DualMap* const & map)
  {
    if (typeDualMap_) delete typeDualMap_;
    typeDualMap_ =  map;
  }

  void setTypeDualMap(DualMap* const & map) { typeDualMap_ =  map; }
  

  const DualMap* getTypeDualMap() const { return typeDualMap_; }


  void replaceStrToPredTemplateMapAndPredDualMap(StrToPredTemplateMap* const& m,
                                                 DualMap* const & predDualMap)
  {
    if (strToPredTemplateMap_)
    {
      StrToPredTemplateMap::iterator it = strToPredTemplateMap_->begin(); 
      for (; it != strToPredTemplateMap_->end(); it++)
        delete (*it).second; //delete PredTemplate*;
      delete strToPredTemplateMap_;
    }
    if (predDualMap_) delete predDualMap_;
    strToPredTemplateMap_ = m;
    predDualMap_ = predDualMap;
  }


  void setStrToPredTemplateMapAndPredDualMap(StrToPredTemplateMap* const & map,
                                             DualMap* const & predDualMap)
  { strToPredTemplateMap_ = map; predDualMap_ = predDualMap; }


  const StrToPredTemplateMap* getStrToPredTemplateMap() const
  { return strToPredTemplateMap_; }


  const DualMap* getPredDualMap() const { return predDualMap_; }


  void replaceStrToFuncTemplateMapAndFuncDualMap(StrToFuncTemplateMap* const& m,
                                                 DualMap* const & funcDualMap)
  {
    if (strToFuncTemplateMap_)
    {
      StrToFuncTemplateMap::iterator it2 = strToFuncTemplateMap_->begin(); 
      for (; it2 != strToFuncTemplateMap_->end(); it2++)
        delete (*it2).second; //delete FuncTemplate*;
      delete strToFuncTemplateMap_;
    }
    if (funcDualMap_) delete funcDualMap_;
    strToFuncTemplateMap_ = m; 
    funcDualMap_ = funcDualMap;
  }


  void setStrToFuncTemplateMapAndFuncDualMap(StrToFuncTemplateMap* const& m,
                                             DualMap* const & funcDualMap)
  { strToFuncTemplateMap_ = m;  funcDualMap_ = funcDualMap; }


  const StrToFuncTemplateMap* getStrToFuncTemplateMap() const
  { return strToFuncTemplateMap_; }


  const DualMap* getFuncDualMap() const { return funcDualMap_; }

  
  void replaceEqualPredTemplate(PredicateTemplate* const & t)
  {
    if (equalPredTemplate_) delete equalPredTemplate_;
    equalPredTemplate_ = t;
  }

  void setEqualPredTemplate(PredicateTemplate* const & t) 
  { equalPredTemplate_ = t; }

  
  const PredicateTemplate* getEqualPredTemplate() const
  { return equalPredTemplate_; }

  void replaceEmptyPredTemplate(PredicateTemplate* const & t)
  {
    if (emptyPredTemplate_) delete emptyPredTemplate_;
    emptyPredTemplate_ = t;
  }

  void setEmptyPredTemplate(PredicateTemplate* const & t) 
  { emptyPredTemplate_ = t; }

  
  const PredicateTemplate* getEmptyPredTemplate() const
  { return emptyPredTemplate_; }

  void replaceEmptyFuncUnaryTemplate(FunctionTemplate* const & t)
  {
    if (emptyFuncUnaryTemplate_) delete emptyFuncUnaryTemplate_;
    emptyFuncUnaryTemplate_ = t;
  }

  void replaceEmptyFuncBinaryTemplate(FunctionTemplate* const & t)
  {
    if (emptyFuncBinaryTemplate_) delete emptyFuncBinaryTemplate_;
    emptyFuncBinaryTemplate_ = t;
  }

  void setEmptyFuncUnaryTemplate(FunctionTemplate* const & t) 
  { emptyFuncUnaryTemplate_ = t; }

  void setEmptyFuncBinaryTemplate(FunctionTemplate* const & t) 
  { emptyFuncBinaryTemplate_ = t; }

  
  const FunctionTemplate* getEmptyFuncUnaryTemplate() const
  { return emptyFuncUnaryTemplate_; }

  const FunctionTemplate* getEmptyFuncBinaryTemplate() const
  { return emptyFuncBinaryTemplate_; }

  void setConstDualMap(ConstDualMap* const & map) { constDualMap_ = map; }

  const ConstDualMap* getConstDualMap() const { return constDualMap_; }

  void setConstantsByType(Array<Array<int>*>* const & cbt) 
  { constantsByType_ = cbt; }

  
  const Array<Array<int>*>* getConstantsByType() const 
  { return constantsByType_; }

  const Array<Array<int>*>* getExternalConstantsByType() const 
  { return externalConstantsByType_; }

  void replaceFuncSet(FunctionSet* const & funcSet)
  {
    if (funcSet_)
    {
      FunctionSet::iterator fit;
      while (!funcSet_->empty())
      { 
        fit = funcSet_->begin();
        funcSet_->erase(fit);
        delete *fit;
      }
      delete funcSet_;
    }
    funcSet_ = funcSet;
  }


  void setFuncSet(FunctionSet* const & funcSet)  { funcSet_ = funcSet; }


  const FunctionSet* getFuncSet() const { return funcSet_; }

  void deleteDB();
  Database* getDB() const { return db_; }
  void setDB(Database* const & db) { db_ = db; }

  void newTrueFalseGroundingsStore();
  TrueFalseGroundingsStore* getTrueFalseGroundingsStore() const
  { return trueFalseGroundingsStore_; }
  void setTrueFalseGroundingsStore(TrueFalseGroundingsStore* const & tfgs)
  { trueFalseGroundingsStore_ = tfgs; }


  void compress();

    // common
  void updatePerOldToNewIds(MLN* const & mln,
                            hash_map<int,int> & oldToNewConstIds);

    // for parser
  void reorderConstants(MLN* const & mln,
                        hash_map<int, PredicateHashArray*>& predIdToPredsMap);

    // for domain0 - postprocess
  void reorderConstants(MLN* const & mln);

    // for domain1-N - postprocess
  void reorderConstants(ConstDualMap* const & map,
                        Array<Array<int>*>* const & cbt,  
                        Array<Array<int>*>* const & ecbt, MLN* const & mln);


    //Caller is responsible for deleting returned pointer
  Predicate* createPredicate(const int& predId, 
                             const bool& includeEqualPreds) const;

    //Caller is responsible for deleting Array and its contents
  void createPredicates(Array<Predicate*>* const & preds,
                        const bool& includeEqualPreds) const;

  void createPredicates(Array<Predicate*>* const & preds,
                        const Array<string>* const & predNames);

  //////////////////////// type functions //////////////////////

    // Returns id of type or -1 if type has been added before.
    // Type id increases by one each time addType() is called.
    // Caller is responsible for deleting name if required. 
  int addType(const char* const& name) 
  { 
    int typeId = typeDualMap_->insert(name); 
    if (typeId < 0) return -1;

    if (typeId != constantsByType_->size())
    {
      cout << "Error: In Domain::addType(). Expected typeId " << typeId 
           << " to be equal to " << constantsByType_->size() << endl;
      exit(-1);
    }
      
    constantsByType_->append(new Array<int>);
    externalConstantsByType_->append(new Array<int>);
    return typeId;
  }


  int getNumTypes() const { return typeDualMap_->getNumInt(); }


    // Caller should NOT delete the returned array.
  const Array<const char*>* getTypeNames() const  
  { return typeDualMap_->getIntToStrArr(); }


    // Caller is responsible for deleting name if required.
    // Returns -1 if type does not exist.
  int getTypeId(const char* const & name) const 
  { return typeDualMap_->getInt(name); }
  
    // Caller should not delete returned const char*
    // Returns NULL if type id does not exist.
  const char* getTypeName(const int typeId) const
  { return typeDualMap_->getStr(typeId); }

  bool isType(const char* const & name) const { return (getTypeId(name) >= 0); }
  bool isType(const int& id) const { return (getTypeName(id) != NULL); }


  //////////////////////// constant functions //////////////////////
    // Caller must not delete returned pointer.
  const Array<int>* getConstantsByType(const int& typeId) const
  {
    return (*constantsByType_)[typeId];
  }

    // Caller must not delete returned pointer.
  const Array<int>* getConstantsByType(const char* const & typeName) const
  {
    int typeId = getTypeId(typeName);
    if (typeId < 0) return NULL;
    return getConstantsByType(typeId);
  }

    // Caller must delete returned pointer.
  const Array<int>* getConstantsByTypeWithExt(const int& typeId) const
  {
    Array<int>* cbt = (*constantsByType_)[typeId];
    Array<int>* retArr = new Array<int>;
    for (int i = 0; i < cbt->size(); i++)
      retArr->append((*cbt)[i]);
    Array<int>* ecbt = (*externalConstantsByType_)[typeId];
    for (int i = 0; i < ecbt->size(); i++)
      retArr->append((*ecbt)[i]);
    retArr->quicksort();
    return retArr;
  }

    // Caller must delete returned pointer.
  const Array<int>* getConstantsByTypeWithExt(const char* const & typeName) const
  { 
    int typeId = getTypeId(typeName);
    if (typeId < 0) return NULL;
    return getConstantsByTypeWithExt(typeId);
  }

	// Returns the index of the given constant of all constants of this type

  int getConstantIndexInType(const int& constId, const int& typeId) const
  {
  	const Array<int>* constArray = getConstantsByType(typeId);
  	for (int i = 0; i < constArray->size(); i++)
  	  if (constId == (*constArray)[i])
        return i;
  	return -1;
  }

  /**
   * Returns the number of constants of a certain type in this domain.
   * External constants are not included.
   * 
   * @param typeId Id of type for which number of constants is retrieved.
   * 
   * @return Number of constants with type typeId in this domain (not including
   * external).
   */  
  int getNumConstantsByType(const int& typeId) const
  { 
    if (!(*constantsByType_)[typeId]) return 0;
    return (*constantsByType_)[typeId]->size();
  }

  /**
   * Returns the number of constants of a certain type in this domain.
   * External constants are not included.
   * 
   * @param typeName Name of type for which number of constants is retrieved.
   * 
   * @return Number of constants of type typeName in this domain (not including
   * external).
   */  
  int getNumConstantsByType(const char* const & typeName) const
  {
    int typeId = getTypeId(typeName);
    if (typeId < 0) return 0;
    return getNumConstantsByType(typeId);
  }

  /**
   * Returns the number of constants of a certain type in this domain.
   * External constants are included.
   * 
   * @param typeId Id of type for which number of constants is retrieved.
   * 
   * @return Number of constants with type typeId in this domain (including external).
   */  
  int getNumConstantsByTypeWithExt(const int& typeId) const
  {
    return ((*constantsByType_)[typeId])->size() +
           ((*externalConstantsByType_)[typeId])->size();
  }


  /**
   * Returns the number of constants of a certain type in this domain.
   * External constants are included.
   * 
   * @param typeName Name of type for which number of constants is retrieved.
   * 
   * @return Number of constants of type typeName in this domain (including external).
   */  
  int getNumConstantsByTypeWithExt(const char* const & typeName) const
  {
    int typeId = getTypeId(typeName);
    if (typeId < 0) return 0;
    return getNumConstantsByTypeWithExt(typeId); 
  }

  /**
   * Adds a constant to the domain.
   * 
   * @param constName Name of constant as a char*
   * @param typeName Type of constant as a char*
   */
  int addConstant(const char* const & constName, const char* const & typeName)
  {
    return addConstant(constName, typeName, false);
  }

  /**
   * Adds an external constant to the domain. External constants are not used when
   * generating groundings of preds.
   * 
   * @param constName Name of constant as a char*
   * @param typeName Type of constant as a char*
   */
  int addExternalConstant(const char* const & constName,
                          const char* const & typeName)
  {
    return addConstant(constName, typeName, true);
  }

    // Returns id of constant or -1 if constant has been added before.
    // Constant id increases by one each time addConstant() is called.
    // Caller is responsible for deleting name and typeName if needed.
  int addConstant(const char* const & constName, const char* const & typeName,
                  const bool& external)
  {
    int typeId =  typeDualMap_->getInt(typeName);
    if (typeId < 0)
    {
      cout << "Warning: failed to add constant " << constName 
           << " because its type of " << typeName << " doesn't exist" << endl;
      return typeId;
    }

      // If constant already exists, then add it to new type
    int constId = getConstantId(constName);
    if (constId >= 0)
    {
      bool newType = true;
      Array<const char*>* prevTypes = getConstantTypeNames(constId);
      for (int i = 0; i < prevTypes->size(); i++)
      {
        if (strcmp((*prevTypes)[i], typeName) == 0)
        {
            // New type
          newType = false;
          break;
        }
      }
      
      if (!newType)
      {
          // Constant already belongs to this type
        return -1;
      }
    }
      // Insert constant into constDualMap_
    constId = constDualMap_->insert(constName, typeId);
    if (constId < 0)
    {
      cout << "Warning: failed to add constant " << constName << " of type " 
           << typeName << endl;
      return constId;
    }

      // Insert id into by type index
    if (external)
    {
      (*externalConstantsByType_)[typeId]->append(constId);
    }
    else
    {
      (*constantsByType_)[typeId]->append(constId);
    }
    return constId;
  }

  /**
   * Replaces the type of a constant from the domain.
   * 
   * @param constName Name of constant being changed.
   * @param oldTypeName Type of constant being chnaged.
   * @param newTypeName New type of the constant.
   * 
   * @return id of constant changed, or -1 if not found.
   */
/*
  int replaceTypeOfConstant(const char* const & constName,
  							const char* const & oldTypeName,
  							const char* const & newTypeName)
  {
    int oldTypeId = typeDualMap_->getInt(oldTypeName);
    if (oldTypeId < 0)
    {
      cout << "Warning: failed to replace type of constant " << constName 
           << " because its type of " << oldTypeName << " doesn't exist" << endl;
      return oldTypeId;
    }
    int newTypeId = typeDualMap_->getInt(newTypeName);
    if (newTypeId < 0)
    {
      cout << "Warning: failed to replace type of constant " << constName 
           << " because its type of " << newTypeName << " doesn't exist" << endl;
      return newTypeId;
    }
    
    int constId = constDualMap_->insert(constName, newTypeId);
    if (constId < 0)
    {
      cout << "Warning: failed to remove constant " << constName << " of type " 
           << oldTypeName << endl;
      return constId;
    }
    
    (*constantsByType_)[oldTypeId]->removeItem(constId);
    (*constantsByType_)[newTypeId]->append(constId);
    return constId;
  }
*/
  
  /**
   * Returns number of (internal) constants.
   */
/*
  int getNumConstants() const
  {
    int count = 0;
    for (int i = 0; i < constantsByType_->size(); i++)
    {
      for (int j = 0; j < (*constantsByType_)[i]->size(); j++)
        count++;
    }
    return count;
  }
*/

  /**
   * Returns number of internal and external constants.
   */
/*
  int getNumConstantsWithExt() const
  {
    return constDualMap_->getNumInt();
  }
*/

  /**
   * Returns the name of a constant in this domain associated with a given
   * id. Caller should not delete returned const char*
   * 
   * @param id Id of constant whose name is to be retrieved.
   * 
   * @return name of constant if id is found as constant or external constant,
   * otherwise NULL.
   */
  const char* getConstantName(const int& id) const
  {
    return constDualMap_->getStr(id);
  }


    // Caller is responsible for deleting name if required.
    // Returns -1 if the constant isn't found
  int getConstantId(const char* const & name) const
  {
    return constDualMap_->getInt(name);
  }


  bool isConstant(const char* const & name) const
  { return (getConstantId(name) >= 0); }
  bool isConstant(const int& id) const { return (getConstantName(id) != NULL); }

    // Caller should delete the constName argument if required
  Array<int>* getConstantTypeIds(const char* const & constName) const 
  { 
    return constDualMap_->getInt2(constName);
  }


  Array<int>* getConstantTypeIds(const int& constId) const
  {
    return constDualMap_->getInt2(constId);
  }


    // Caller should not delete returned const char*
    // Returns NULL if id does not exist.
  Array<const char*>* getConstantTypeNames(const int& constId) const
  { 
    Array<int>* typeIds = getConstantTypeIds(constId);
    if (typeIds == NULL) return NULL;
    else
    {
      Array<const char*>* typeNames = new Array<const char*>;
      for (int i = 0; i < typeIds->size(); i++)
        typeNames->append(getTypeName((*typeIds)[i]));
      return typeNames;
    }
  }




  //////////////////////// predicate functions ///////////////////

    // Returns id of predicate or -1 if PredicateTemplate has been added before.
    // predTemplate should have its name_, termTypes_, and domain_
    // assigned before this function is called.
    // Predicate id increases by one each time addPredicateTemplate() is called.
    // Caller should not delete predTemplate. 
  int addPredicateTemplate(const PredicateTemplate* const & predTemplate)
  {
    int predId = predDualMap_->insert(predTemplate->getName());
    if (predId < 0)
    {
      cout << "Warning: failed to add predicate template " 
           << predTemplate->getName() << endl;
      return predId;
    }

    const char* predName = predDualMap_->getStr(predId);
     //strToPredTemplateMap_ shares the predName object with predDualMap_
    (*strToPredTemplateMap_)[predName] = predTemplate; 
    return predId;
  }


  int getNumPredicates() const { return predDualMap_->getNumInt(); }


    // Caller should NOT delete the returned array and its contents.
  const Array<const char*>* getPredicateNames() const 
  { return predDualMap_->getIntToStrArr(); }


  void getNonEqualPredicateNames(Array<string>& predNames) const
  {
    for (int i = 0; i < getNumPredicates(); i++)
    {
      if (getPredicateTemplate(i)->isEqualPredicateTemplate()) continue;
      if (getPredicateTemplate(i)->isInternalPredicateTemplate()) continue;
      if (getPredicateTemplate(i)->isPredicateTemplateFromFunction()) continue;
      predNames.append(getPredicateName(i));
    }
  }


    // Caller is responsible for deleting name if required.
    // Returns id of predicate if it exist, otherwise returns -1
  int getPredicateId(const char* const & name) const
  { return predDualMap_->getInt(name); }

    // Caller should not delete returned const char*
    // Returns name of predicate if it exists, otherwise returns NULL
  const char* getPredicateName(const int& id) const
  { return predDualMap_->getStr(id); }

    //Returns the term type ids of the predicate with the specified name, or
    //NULL if there is no predicate with specified name.
    //Caller should delete name if required.
    //Caller should not delete the returned const Array<int>*
  const Array<int>* getPredicateTermTypesAsInt(const char* const & name) const
  {
    StrToPredTemplateMap::iterator it;
    if ((it=strToPredTemplateMap_->find(name)) == strToPredTemplateMap_->end())
      return NULL;
    return (*it).second->getTermTypesAsInt();
  }


    //Returns the term type ids of the predicate with the specified id, or
    //NULL if there is no predicate with specified id.
    //Caller should not delete the returned const Array<int>*
  const Array<int>* getPredicateTermTypesAsInt(const int& id) const
  { return getPredicateTermTypesAsInt(getPredicateName(id)); }


    //Returns the term type ids of the predicate with the specified name, or
    //NULL if there is no predicate with specified name.
    //Caller should delete name if required.
    //Caller should not delete the returned const Array<const char*>*
  const Array<const char*>* getPredicateTermTypesAsStr(const char* const& name)
    const
  {
    StrToPredTemplateMap::iterator it;
    if ((it=strToPredTemplateMap_->find(name)) == strToPredTemplateMap_->end())
      return NULL;
    return (*it).second->getTermTypesAsStr();
  }


    //Returns the term type ids of the predicate with the specified id, or
    //NULL if there is no predicate with specified id.
    //Caller should not delete the returned const Array<const char*>*
  const Array<const char*>* getPredicateTermTypesAsStr(const int& id)
  { return getPredicateTermTypesAsStr(getPredicateName(id)); }


    // Returns the PredicateTemplate* with the given name or NULL if there is no
    // PredicateTemplate with the specified name.
    // Caller is responsible for deleting name if required.
    // Caller should not delete the returned PredicateTemplate* nor modify it.
  const PredicateTemplate* getPredicateTemplate(const char* const & name) const
  {
    StrToPredTemplateMap::iterator it;
    if ((it=strToPredTemplateMap_->find(name)) == strToPredTemplateMap_->end())
      return NULL;
    return (*it).second;
  }


    // Returns the PredicateTemplate* with the given id.
    // Caller should not delete the returned PredicateTemplate* nor modify it.
  const PredicateTemplate* getPredicateTemplate(const int& id) const
  {
    const char* predName = ((Domain*)this)->getPredicateName(id);
    return getPredicateTemplate(predName);
  }

 /**
  * Finds the maximum arity of all predicates present in the domain.
  *  
  * @return maximum arity of all predicates in the domain.
  */
  const int getHighestPredicateArity() const
  {
  	int highestArity = 1;
  	int arity;
  	for (int i = 0; i < getNumPredicates(); i++)
    {
      arity = getPredicateTemplate(i)->getNumTerms();
      if (arity > highestArity) highestArity = arity;
    }
    
  	return highestArity;
  }

  const PredicateTemplate* getEqualPredicateTemplate() const
  { return equalPredTemplate_; }
  
  const PredicateTemplate* getEmptyPredicateTemplate() const
  { return emptyPredTemplate_; }
 
  const FunctionTemplate* getEmptyFunctionUnaryTemplate() const
  { return emptyFuncUnaryTemplate_; }

  const FunctionTemplate* getEmptyFunctionBinaryTemplate() const
  { return emptyFuncBinaryTemplate_; }

  bool isPredicate(const char* const & name) const 
  {
  	return (getPredicateId(name)>=0 ||
  			strcmp(name, PredicateTemplate::GT_NAME) == 0 || 
  			strcmp(name, PredicateTemplate::LT_NAME) == 0 ||
  			strcmp(name, PredicateTemplate::GTEQ_NAME) == 0 ||
  			strcmp(name, PredicateTemplate::LTEQ_NAME) == 0 ||
  			strcmp(name, PredicateTemplate::SUBSTR_NAME) == 0);
  	
  }

  bool isPredicate(const int& id) const {return (getPredicateName(id) != NULL);}


  //////////////////////// function functions ///////////////////

    // Returns id of function or -1 if function has been added before.
    // functionTemplate should have its name_, termTypes_, domain_, 
    // and retTypeName_ assigned before this function is called.
    // Function id increases by one each time addFunctionTemplate() is called.
    // Caller should not delete funcTemplate. 
  int addFunctionTemplate(const FunctionTemplate* const & funcTemplate)
  {
    int funcId = funcDualMap_->insert(funcTemplate->getName());
    if (funcId < 0)
    {
      cout << "Warning: failed to add function template " 
           << funcTemplate->getName() << endl;
      return funcId;
    }

    const char* funcName = funcDualMap_->getStr(funcId);
     //strToFuncTemplateMap_ shares the funcName object with funcDualMap_
    (*strToFuncTemplateMap_)[funcName] = funcTemplate; 
    return funcId;
  }


  int getNumFunctions() const { return funcDualMap_->getNumInt(); }


    // Caller should NOT delete the returned array but not its contents.
  const Array<const char*>* getFunctionNames() const 
  { return funcDualMap_->getIntToStrArr(); }


    // Caller is responsible for deleting name if required.
  int getFunctionId(const char* const & name) const
  { return funcDualMap_->getInt(name); }

    // Caller should not delete returned const char*
  const char* getFunctionName(const int& id) const
  { return funcDualMap_->getStr(id); }

    //Returns the term type ids of the function with the specified name, or
    //NULL if there is no function with specified name.
    //Caller should delete name if required.
    //Caller should not delete the returned const Array<int>*
  const Array<int>* getFunctionTermTypesAsInt(const char* const & name) const
  {
    StrToFuncTemplateMap::iterator it;
    if ((it=strToFuncTemplateMap_->find(name)) == strToFuncTemplateMap_->end())
      return NULL;
    return (*it).second->getTermTypesAsInt();
  }


    //Returns the term type ids of the function with the specified name, or
    //NULL if there is no function with specified name.
    //Caller should delete name if required.
    //Caller should not delete the returned const Array<const char*>*
  const Array<const char*>* getFunctionTermTypesAsStr(const char* const & name)
    const
  {
    StrToFuncTemplateMap::iterator it;
    if ((it=strToFuncTemplateMap_->find(name)) == strToFuncTemplateMap_->end())
      return NULL;
    return (*it).second->getTermTypesAsStr();
  }


    // Returns the FunctionTemplate with the given name 
    // Caller is responsible for deleting name if required.
    // Caller should not delete the returned FunctionTemplate* nor modify it.
  const FunctionTemplate* getFunctionTemplate(const char* const & name) const
  {
    StrToFuncTemplateMap::iterator it;
    if ((it=strToFuncTemplateMap_->find(name)) == strToFuncTemplateMap_->end())
      return NULL;
    return (*it).second;
  }

    // Returns the FunctionTemplate* with the given id or NULL if there is no
    // FunctionTemplate with the specified id.
    // Caller should not delete the returned FunctionTemplate* nor modify it.
  const FunctionTemplate* getFunctionTemplate(const int& id) const
  {
    const char* funcName = ((Domain*)this)->getFunctionName(id);
    return getFunctionTemplate(funcName);
  }
  
  
  bool isFunction(const char* const & name) const 
  {
  	return (getFunctionId(name)>=0 ||
  			strcmp(name, FunctionTemplate::SUCC_NAME) == 0 || 
  			strcmp(name, FunctionTemplate::PLUS_NAME) == 0 ||
  			strcmp(name, FunctionTemplate::MINUS_NAME) == 0 ||
  			strcmp(name, FunctionTemplate::TIMES_NAME) == 0 ||
  			strcmp(name, FunctionTemplate::DIVIDEDBY_NAME) == 0 ||
  			strcmp(name, FunctionTemplate::MOD_NAME) == 0 ||
  			strcmp(name, FunctionTemplate::CONCAT_NAME) == 0);
  }
  
  bool isFunction(const int& id) const { return (getFunctionName(id) != NULL); }


    // Caller shouldn't delete the returned function set nor modify its contents
  const FunctionSet* getFunctionMappings() const { return funcSet_; }

  
    // f is owned by Domain and caller should not delete it
    // returns true if f is added, and false if f already exists
  bool addFunctionMapping(Function* f) 
  {
  	assert(f->getRetConstId() >= 0);
    if (funcSet_->find(f) ==  funcSet_->end())
    {
      funcSet_->insert(f);
      return true;
    }
    return false;
  }

  
  int getFunctionRetConstId(Function* const & f) const
  {
    FunctionSet::iterator it;
    if ((it = funcSet_->find(f)) == funcSet_->end()) return -1;
    return (*it)->getRetConstId();
  }


  void printPredicateTemplates(ostream& out) const
  {
    StrToPredTemplateMap::iterator it = strToPredTemplateMap_->begin();
    for (; it !=  strToPredTemplateMap_->end(); it++)
    {
//      if (!(*it).second->isEqualPredicateTemplate() &&
//      	  !(*it).second->isInternalPredicateTemplate() &&
//      	  !(*it).second->isPredicateTemplateFromFunction())
      if (!(*it).second->isEqualPredicateTemplate() &&
      	  !(*it).second->isInternalPredicateTemplate())
        out << *((*it).second) << endl;
    }
  }

  /**
   * print out the function declarations.
   * 
   * @out output stream to which the declarations are written
   */
  void printFunctionTemplates(ostream& out) const
  {
    StrToFuncTemplateMap::iterator it = strToFuncTemplateMap_->begin();
    for (; it !=  strToFuncTemplateMap_->end(); it++)
    {
      if (!(*it).second->isInternalFunctionTemplate())
      	out << *((*it).second) << endl;
    }
  }

  int getNumNonEvidenceAtoms() const;
  
  Predicate* getNonEvidenceAtom(const int& index) const;
  
  int addPredBlock(Predicate* const & predBlock) const;

  int getBlock(const Predicate* const & pred) const;
  
  int getBlock(const GroundPredicate* const & pred) const;
  
  int getNumPredBlocks() const;
  
  const Array<Predicate*>* getPredBlocks() const;

  const Predicate* getPredBlock(const int& index) const;
  
  /**
   * Caller should not delete returned Predicate.
   * 
   * @param index Index of the block.
   */
  const Predicate* getTruePredInBlock(const int& index) const
  {
    assert(index <= truePredsInBlock_->size());
    return (*truePredsInBlock_)[index];
  }

  /**
   * Sets the true predicate in a block. The truth value must also be set in
   * the database.
   * 
   * @param index Index of the block.
   * @param pred Predicate being set to the true one in the block.
   */
  void setTruePredInBlock(const int& index, Predicate* const & pred) const
  {
    assert(index <= truePredsInBlock_->size());
    if ((*truePredsInBlock_)[index]) delete (*truePredsInBlock_)[index];
    (*truePredsInBlock_)[index] = pred;
  }

  const int getBlockSize(const int& index) const
  {
    return (*blockSizes_)[index];
  }

  const Array<bool>* getBlockEvidenceArray() const;
  
  const bool getBlockEvidence(const int& index) const;

  void setBlockEvidence(const int& index, const bool& value) const;

  /**
   * Picks a predicate at random from a block. Caller should delete the
   * returned Predicate*.
   * 
   * @param index Index of block from which to generate the Predicate.
   */
  const Predicate* getRandomPredInBlock(const int& index) const
  {
    const int chosen = random() % getBlockSize(index);
    return getPredInBlock(chosen, index);
  }
  
  /**
   * Gets one ground predicate from a block. Caller should delete the
   * returned Predicate*.
   * 
   * @param grounding Index of grounding in block to generate.
   * @param index Index of block from which to generate the Predicate.
   */
  const Predicate* getPredInBlock(const int& grounding, const int& index) const
  {
    Predicate* pred = (Predicate*)getPredBlock(index);
    assert(grounding <= getBlockSize(index));
    Array<Predicate*> predArr;
    pred->getGroundingNumber(this, predArr, grounding);
    assert(predArr.size() == 1);
    const Predicate* gndPred = new Predicate(*predArr[0]);
    predArr.deleteItemsAndClear();
    return gndPred;
  }
  
  /**
   * Gets the index in a block of a predicate
   * 
   * @param pred Predicate being searched for in the block.
   * @param blockIdx Index of block to search.
   * 
   * @return Index in the block of the predicate, if present; otherwise, -1.
   */
  const int getIndexOfPredInBlock(const Predicate* const & pred,
                                  const int& blockIdx) const
  {
    assert(((Predicate*)pred)->isGrounded());

    for (int i = 0; i < getBlockSize(blockIdx); i++)
    {
      const Predicate* predInBlock = getPredInBlock(i, blockIdx);
      if (((Predicate*)pred)->same((Predicate*)predInBlock))
      {
        delete predInBlock;
        return i;
      }
      delete predInBlock;
    }
    return -1;
  }

  /**
   * Computes the number of non-evidence atoms for each first-order predicate.
   */
  void computeNumNonEvidAtoms();
   
  void setNumTrueNonEvidGndingsPerClause(Array<double>*
                                         numTrueNonEvidGndingsPerClause)
  {
    numTrueNonEvidGndingsPerClause_ = numTrueNonEvidGndingsPerClause;
  }
  
  void setNumFalseNonEvidGndingsPerClause(Array<double>*
                                          numFalseNonEvidGndingsPerClause)
  {
    numFalseNonEvidGndingsPerClause_ = numFalseNonEvidGndingsPerClause;
  }

  void setNumTrueNonEvidGndings(const int & clauseIdx, const double & count)
  {
    if (numTrueNonEvidGndingsPerClause_->size() <= clauseIdx)
      numTrueNonEvidGndingsPerClause_->growToSize(clauseIdx + 1);
    (*numTrueNonEvidGndingsPerClause_)[clauseIdx] = count;
  }

  void setNumFalseNonEvidGndings(const int & clauseIdx, const double & count)
  {
    if (numFalseNonEvidGndingsPerClause_->size() <= clauseIdx)
      numFalseNonEvidGndingsPerClause_->growToSize(clauseIdx + 1);
    (*numFalseNonEvidGndingsPerClause_)[clauseIdx] = count;
  }

  const double getNumTrueNonEvidGroundings(const int & clauseIdx) const
  {
    return (*numTrueNonEvidGndingsPerClause_)[clauseIdx];
  }  

  const double getNumFalseNonEvidGroundings(const int & clauseIdx) const
  {
    return (*numFalseNonEvidGndingsPerClause_)[clauseIdx];
  }  

  const double getNumNonEvidGroundings(const int & clauseIdx) const
  {
    return ((*numTrueNonEvidGndingsPerClause_)[clauseIdx] +
            (*numFalseNonEvidGndingsPerClause_)[clauseIdx]);
  }  

    //get the predicate corresponding to the tuple of constants
    //Note: caller is responsible for deleting it
  Predicate * getPredicate(Array<int> * const & constants, int predId)
  {
    Term *term;
    const PredicateTemplate *pt = getPredicateTemplate(predId);
    Predicate *pred = new Predicate(pt);
    for (int i = 0; i < constants->size(); i++)
    {
      term = new Term((*constants)[i]);
      pred->appendTerm(term);
    }
    return pred;
  }

    //get the predicate corresponding to the predId
    //Note: caller is responsible for deleting it
  Predicate * getPredicate(const int& predId)
  {
    const PredicateTemplate *pt = getPredicateTemplate(predId);
    Predicate *pred = new Predicate(pt);
    return pred;
  }

 private:
 
  void changePredTermsToNewIds(Predicate* const & p,
                               hash_map<int,int>& oldToNewConstIds);


 private:
  DualMap* typeDualMap_;  //maps type id to name and vice versa
  ConstDualMap* constDualMap_; //maps constant id to name and vice versa
  DualMap* predDualMap_;  //maps predicate id to name and vice versa
  DualMap* funcDualMap_;  //maps func id to name and vice versa
    //maps predicate name to PredicateTemplate*. shares key with predDualMap_
  StrToPredTemplateMap* strToPredTemplateMap_;
    //maps function name to FunctionTemplate*. shares its key with funcDualMap_
  StrToFuncTemplateMap* strToFuncTemplateMap_;

    // corresponds to the predicate EQUAL_NAME(SK_ANY, SK_ANY)
  PredicateTemplate* equalPredTemplate_;
  
    // empty binary predicate template used when infix operators are used
  PredicateTemplate* emptyPredTemplate_;
  
    // empty unary function template used when infix operators are used
  FunctionTemplate* emptyFuncUnaryTemplate_;

    // empty binary function template used when infix operators are used
  FunctionTemplate* emptyFuncBinaryTemplate_;
  
    // constantsByType_[0] is an Array<int> of constants of type 0
  Array<Array<int>*>* constantsByType_;
    // externalConstantsByType_[0] is an Array<int> of ext. constants of type 0
  Array<Array<int>*>* externalConstantsByType_;
  Database* db_;
  TrueFalseGroundingsStore* trueFalseGroundingsStore_; //for sampling clauses
  FunctionSet* funcSet_;
    // (Partially grounded) predicates which represent one block. Terms which
    // are still variables are the vars with !
  Array<Predicate*>* predBlocks_;
    // Ground predicate in each block which is true
  Array<Predicate*>* truePredsInBlock_;
    // Grounding index in each block which is true (i-th combination of
    // constants)
  Array<int>* trueGroundingsInBlock_;
    // Size of each block
  Array<int>* blockSizes_;
    // Flags indicating if block is fulfilled by evidence
  Array<bool>* blockEvidence_;
  
    // numNonEvidAtomsPerPred_[p] holds the number of non-evidence groundings
    // of the first-order predicate with index p
  Array<int>* numNonEvidAtomsPerPred_;

    // numTrueNonEvidGndingsPerClause_[c] holds the total number of true
    // groundings of the first-order clause not satisfied by evidence
  Array<double>* numTrueNonEvidGndingsPerClause_;
  Array<double>* numFalseNonEvidGndingsPerClause_;
};


#endif
