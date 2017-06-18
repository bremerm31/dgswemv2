#ifndef CLASS_PROBLEM_H
#define CLASS_PROBLEM_H

#include "../../class_mesh.h"
#include "../../general_definitions.h"

class PROBLEM {
private:
    MESH* mesh;
    
	double t = 0;
    double dt;

public:
    PROBLEM();
    ~PROBLEM();

	void Solve(int, double, double, double);

private:
	void WriteDataVTK();

	void EETimeStepper(int);
	void RK2TimeStepper(int);
	void RK3TimeStepper(int);
	void RK4TimeStepper(int);

	void Timestep();

    void ComputeUVA(Element<>*);
    void ComputeF(Element<>*);
	void ComputeS(Element<>*);

	void InterfaceFlowAverage(INTERFACE*);
	void ComputeBoundaryInterfaceF(INTERFACE*);
	void LLFNumericalFlux(INTERFACE*);

    void OceanInterfaceSetBC(INTERFACE*);	
    void LandInterfaceSetBC(INTERFACE*);
    void FlowInterfaceSetBC(INTERFACE*);
//	void RadiationInterfaceSetBC(INTERFACE*);
//  void NodeInterfaceReceiveData(INTERFACE*);
//  void NodeInterfaceSendData(INTERFACE*);

};

#endif