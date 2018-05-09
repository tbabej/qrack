//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano 2017, 2018. All rights reserved.
//
// This is a multithreaded, universal quantum register simulation, allowing
// (nonphysical) register cloning and direct measurement of probability and
// phase, to leverage what advantages classical emulation of qubits can have.
//
// Licensed under the GNU General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/gpl-3.0.en.html
// for details.

#include "qengine_cpu.hpp"

#include "common/par_for.hpp"

namespace Qrack {

/// Measurement gate
bool QEngineCPU::M(bitLenInt qubit)
{
    if (runningNorm != 1.0) {
        NormalizeState();
    }

    bool result;
    double prob = Rand();
    double angle = Rand() * 2.0 * M_PI;
    double cosine = cos(angle);
    double sine = sin(angle);
    Complex16 nrm;

    bitCapInt qPowers = 1 << qubit;
    double oneChance = Prob(qubit);

    result = (prob < oneChance) && oneChance > 0.0;
    double nrmlzr = 1.0;
    if (result) {
        if (oneChance > 0.0) {
            nrmlzr = oneChance;
        }

        nrm = Complex16(cosine, sine) / nrmlzr;

        PAR_FOR(0, maxQPower) {
            if ((idx & qPowers) == 0) {
                stateVec[idx] = Complex16(0.0, 0.0);
            } else {
                stateVec[idx] = nrm * stateVec[idx];
            }
        } PAR_FOR_END;
    } else {
        if (oneChance < 1.0) {
            nrmlzr = sqrt(1.0 - oneChance);
        }

        nrm = Complex16(cosine, sine) / nrmlzr;

        PAR_FOR(0, maxQPower) {
            if ((idx & qPowers) == 0) {
                stateVec[idx] = nrm * stateVec[idx];
            } else {
                stateVec[idx] = Complex16(0.0, 0.0);
            }
        } PAR_FOR_END;
    }

    UpdateRunningNorm();

    return result;
}

// Apply X ("not") gate to each bit in "length," starting from bit index
// "start"
void QEngineCPU::X(bitLenInt start, bitLenInt length)
{
    // First, single bit operations are better optimized for this special case:
    if (length == 1) {
        X(start);
        return;
    }

    // As a fundamental gate, the register-wise X could proceed like so:
    // for (bitLenInt idx = 0; idx < length; idx++) {
    //    X(start + idx);
    //}

    // Basically ALL register-wise gates proceed by essentially the same
    // algorithm as this simple X gate.

    // We first form bit masks for those qubits involved in the operation, and
    // those not involved in the operation. We might have more than one
    // register involved in the operation in general, but we only have one, in
    // this case.
    bitCapInt inOutMask = ((1 << length) - 1) << start;
    bitCapInt otherMask = ((1 << qubitCount) - 1) ^ inOutMask;

    // Sometimes we transform the state in place. Alternatively, we often
    // allocate a new permutation state vector to transfer old probabilities
    // and phases into.
    Complex16 *nStateVec = AllocStateVec(maxQPower);

    // This function call is a parallel "for" loop. We have several variants of
    // the parallel for loop. Some skip certain permutations in order to
    // optimize. Some take a new permutation state vector for output, and some
    // just transform the permutation state vector in place.
    PAR_FOR(0, maxQPower) {
        // Set nStateVec, indexed by the loop control variable (idx) with
        // the X'ed bits inverted, with the value of stateVec indexed by
        // idx.

        // This is the body of the parallel "for" loop. We iterate over
        // permutations of bits.  We're going to transform from input
        // permutation state to output permutation state, and transfer the
        // probability and phase of the input permutation to the output
        // permutation.  These are the bits that aren't involved in the
        // operation.

        // bitCapInt otherRes = (idx & otherMask);

        // These are the bits in the register that is being operated on. In
        // all permutation states, the bits acted on by the gate should be
        // transformed in the logically appropriate way from input
        // permutation to output permutation. Since this is an X gate, we
        // take the involved bits and bitwise NOT them.

        // bitCapInt inOutRes = ((~idx) & inOutMask);

        // Now, we just transfer the untransformed input state's phase and
        // probability to the transformed output state.

        // nStateVec[inOutRes | otherRes] = stateVec[idx];

        // (We can do this all in one line:)
        nStateVec[(idx & otherMask) | ((~idx) & inOutMask)] = stateVec[idx];

        // For other operations, like the quantum equivalent of a logical
        // "AND," we might have two input registers and one output
        // register. The transformation would be that we use bit masks to
        // bitwise "AND" the input values in every permutation and place
        // this logical result into the output register with another bit
        // mask, for every possible permutation state. Basically all the
        // register-wise operations in Qrack proceed this same way.
    } PAR_FOR_END;
    // We replace our old permutation state vector with the new one we just
    // filled, at the end.
    ResetStateVec(nStateVec);
}

/// Bitwise CNOT
void QEngineCPU::CNOT(bitLenInt start1, bitLenInt start2, bitLenInt length)
{
    bitCapInt reg1Mask = ((1 << length) - 1) << start1;
    bitCapInt reg2Mask = ((1 << length) - 1) << start2;
    bitCapInt otherMask = maxQPower - 1;
    otherMask ^= reg1Mask | reg2Mask;

    Complex16 *nStateVec = AllocStateVec(maxQPower);

    PAR_FOR(0, maxQPower) {
        bitCapInt otherRes = (idx & otherMask);
        bitCapInt reg1Res = idx & reg1Mask;
        bitCapInt reg2Res = ((reg1Res >> start1) ^ ((idx & reg2Mask) >> start2)) << start2;
        nStateVec[reg1Res | reg2Res | otherRes] = stateVec[idx];
    } PAR_FOR_END;
    // We replace our old permutation state vector with the new one we just filled, at the end.
    ResetStateVec(nStateVec);
}

/// Bitwise "Anti-"CNOT - NOT operation if control is 0
void QEngineCPU::AntiCNOT(bitLenInt start1, bitLenInt start2, bitLenInt length)
{
    bitCapInt reg1Mask = ((1 << length) - 1) << start1;
    bitCapInt reg2Mask = ((1 << length) - 1) << start2;
    bitCapInt otherMask = maxQPower - 1;
    otherMask ^= reg1Mask | reg2Mask;

    Complex16 *nStateVec = AllocStateVec(maxQPower);

    PAR_FOR(0, maxQPower) {
        bitCapInt otherRes = (idx & otherMask);
        bitCapInt reg1Res = idx & reg1Mask;
        bitCapInt reg2Res = ((((~reg1Res) & reg1Mask) >> start1) ^ ((idx & reg2Mask) >> start2)) << start2;
        nStateVec[reg1Res | reg2Res | otherRes] = stateVec[idx];
    } PAR_FOR_END;
    // We replace our old permutation state vector with the new one we just filled, at the end.
    ResetStateVec(nStateVec);
}

/// Bitwise CCNOT
void QEngineCPU::CCNOT(bitLenInt control1, bitLenInt control2, bitLenInt target, bitLenInt length)
{
    bitCapInt reg1Mask = ((1 << length) - 1) << control1;
    bitCapInt reg2Mask = ((1 << length) - 1) << control2;
    bitCapInt reg3Mask = ((1 << length) - 1) << target;
    bitCapInt otherMask = maxQPower - 1;
    otherMask ^= reg1Mask | reg2Mask | reg3Mask;

    Complex16 *nStateVec = AllocStateVec(maxQPower);

    PAR_FOR(0, maxQPower) {
        bitCapInt otherRes = (idx & otherMask);
        bitCapInt reg1Res = idx & reg1Mask;
        bitCapInt reg2Res = idx & reg2Mask;
        bitCapInt reg3Res = (((reg1Res >> control1) & (reg2Res >> control2)) ^ ((idx & reg3Mask) >> target)) << target;
        nStateVec[reg1Res | reg2Res | reg3Res | otherRes] = stateVec[idx];
    } PAR_FOR_END;
    // We replace our old permutation state vector with the new one we just filled, at the end.
    ResetStateVec(nStateVec);
}

/// Bitwise "Anti-"CCNOT - NOT operation if both control bits are 0
void QEngineCPU::AntiCCNOT(bitLenInt control1, bitLenInt control2, bitLenInt target, bitLenInt length)
{
    bitCapInt reg1Mask = ((1 << length) - 1) << control1;
    bitCapInt reg2Mask = ((1 << length) - 1) << control2;
    bitCapInt reg3Mask = ((1 << length) - 1) << target;
    bitCapInt otherMask = maxQPower - 1;
    otherMask ^= reg1Mask | reg2Mask | reg3Mask;

    Complex16* nStateVec = AllocStateVec(maxQPower);

    PAR_FOR(0, maxQPower) {
        bitCapInt otherRes = (idx & otherMask);
        bitCapInt reg1Res = idx & reg1Mask;
        bitCapInt reg2Res = idx & reg2Mask;
        bitCapInt reg3Res = (((((~reg1Res) & reg1Mask) >> control1) & (((~reg2Res) & reg2Mask) >> control2)) ^ ((idx & reg3Mask) >> target)) << target;
        nStateVec[reg1Res | reg2Res | reg3Res | otherRes] = stateVec[idx];
    } PAR_FOR_END;
    // We replace our old permutation state vector with the new one we just filled, at the end.
    ResetStateVec(nStateVec);
}

/// Bitwise swap
void QEngineCPU::Swap(bitLenInt start1, bitLenInt start2, bitLenInt length)
{
    int distance = start1 - start2;
    if (distance == 0) {
        return;
    }

    bitCapInt reg1Mask = ((1 << length) - 1) << start1;
    bitCapInt reg2Mask = ((1 << length) - 1) << start2;
    bitCapInt otherMask = maxQPower - 1;
    otherMask ^= reg1Mask | reg2Mask;

    Complex16* nStateVec = AllocStateVec(maxQPower);

    PAR_FOR(0, maxQPower) {
        bitCapInt otherRes = (idx & otherMask);
        bitCapInt reg1Res = ((idx & reg1Mask) >> start1) << start2;
        bitCapInt reg2Res = ((idx & reg2Mask) >> start2) << start1;
        nStateVec[reg1Res | reg2Res | otherRes] = stateVec[idx];
    } PAR_FOR_END;
    // We replace our old permutation state vector with the new one we just filled, at the end.
    ResetStateVec(nStateVec);
}

/// Phase flip always - equivalent to Z X Z X on any bit in the QEngineCPU
void QEngineCPU::PhaseFlip()
{
    PAR_FOR(0, maxQPower) {
        stateVec[idx] = -stateVec[idx];
    } PAR_FOR_END;
}

}
