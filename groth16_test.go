package playsnark

import "testing"

func TestGroth16(t *testing.T) {
	r1cs := createR1CS()
	//s := createWitness(r1cs)
	qap := ToQAP(r1cs)
	//diff := qap.nbVars - qap.nbIO
	_ = NewGroth16TrustedSetup(qap)

}
