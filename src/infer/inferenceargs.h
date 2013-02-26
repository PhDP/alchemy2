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
#ifndef INFERENCEARGS_H_
#define INFERENCEARGS_H_

#include "inferutil.h"

/**
 * This file contains the command line arguments for performing inference and
 * the data structures to store them. The default values are declared here with
 * the variable declaration.
 *
 * Naming convention is: The variable should start with an a (for argument),
 * then the abbreviation for the inference algorithm to which the argument is
 * specific (if any), then the variable name.
 */

bool aisQueryEvidence = false;

  // Holds the names of the input mln files comma-separated
char* ainMLNFiles     = NULL;

  // Holds closed-world non-evidence preds
char* aClosedWorldPredsStr = NULL;
  // Holds open-world evidence preds
char* aOpenWorldPredsStr = NULL;

  // MAP Inference
bool  amapPos = false;
bool  amapAll = false;
  // Gibbs Sampling
bool  agibbsInfer = false;
  // MC-SAT
bool  amcsatInfer = false;
  // Simulated Tempering
bool  asimtpInfer = false;
  // Belief propagation
bool  abpInfer = false;
  // Expanding Frontier Belief Propagation
bool  aefbpInfer = false;
  // Decision network: output total utility and assignment of action predicates
bool  adecisionInfer = false;
  // No inference, just output network
bool  aoutputNetwork = false;

  // Lazy state or not?
bool aLazy = false;

  // Hybrid or not?
bool aHybrid = false;

bool aSA = false;

char* aProposalStdev = (char*)"tmp.sd";
bool bOptimum = false;

char*  aContSamples = NULL;

//inference with maxwalksat?
bool aMaxOrNot = false;
  // Seed for inference algorithm
int aSeed = -1;
  // Limit in kbytes before using lazy version
int aMemLimit = -1;
  // If true, atoms won't be deactivated to save space
bool aLazyNoApprox = false;

  // MaxWalkSat params
//int  amwsMaxSteps     = 100000;
int  amwsMaxSteps     = 1000000;
int  amwsTries        = 1;
int  amwsNumSolutions = 10;
int  amwsTargetWt     = 0;
bool amwsHard         = false;
int  amwsHeuristic    = 2; // 0 = RANDOM, 1 = BEST, 2 = TABU, 3 = SampleSat
int  amwsTabuLength   = 10;
bool amwsLazyLowState = false;

  // SampleSat params
int  assSaRatio = 50;       // percent of SA steps
int  assSaTemp  = 80;       // temperature/100: SA temperature
bool assLateSa  = true;    // sa only at a plateur
double aSATempDownRatio = 0.9;

  // MCMC params
int amcmcNumChains    = 10;
int amcmcBurnMinSteps = 100;
int amcmcBurnMaxSteps = 100;
int amcmcMinSteps     = -1;
int amcmcMaxSteps     = 1000;
int amcmcMaxSeconds   = -1;

  // Gibbs params
double agibbsDelta           = 0.05;
double agibbsEpsilonError    = 0.01;
double agibbsFracConverged   = 0.95;
int    agibbsWalksatType     = 1;
bool   agibbsTestConvergence = false;
int    agibbsSamplesPerTest  = 100;

  // Simulated Tempering params
int asimtpSubInterval = 2;
int asimtpNumST       = 3;
int asimtpNumSwap     = 10;

  // Belief Propagation params
bool aliftedInfer = false;
bool auseHC = false;
bool auseCT = false;
double abpConvergenceThresh = 1e-4;
int abpConvergeRequiredItrCnt = 20;
bool aexplicitRep = false;
int ahcCreateType = Basic;
double ahcCreateNoise = 0.0;
int alncIter = 0;
//predicates for which hypercubes should not be created
char *anoHCPredsStr = NULL;

  // Produce clause counts
bool aclauseCounts = false;

#endif /*INFERENCEARGS_H_*/
