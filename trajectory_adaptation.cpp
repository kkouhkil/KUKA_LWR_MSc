#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include
#include
#include
#include
#include
#include
#include
#include
#include
#include
<time.h>
<math.h>
<signal.h>
<iostream>
<Eigen/Dense>
<cstdlib>
<stdlib.h>
<stdio.h>
<unistd.h>
<algorithm>
#include "friComDrv.h"
#include "kukaLWR.h"
#include "model.h"
#define N_MAX 5
#define Ff 0.02
#define TORQUELIMIT 2
using namespace std;
using namespace Eigen;
bool pinv_damped(const MatrixXd A, MatrixXd &invA, VectorXd &sigma, double eps);
typedef enum
{
POSITION,
TORQUE
}Modes;
Modes kuka_mode = TORQUE;
int main(int argc, char *argv[])
{
double t = 0;
double tminus1 = 0;
kukaLWR kuka;
//************************************************************ JOINTS ESSENTIAL VARAIBLES ******************************************************************

float jointPosRadKuka[LBR_MNJ];
float jointPrePosRadKuka[LBR_MNJ];
float jointDesPosRadKuka[LBR_MNJ];
float jointCurDotPosRadKuka[LBR_MNJ];
float jointDesTrq[LBR_MNJ];
float jointLimitedExtTrqEst[LBR_MNJ];

//************************************************************ PROGRAM TIME CYCLE **************************************************************************

float Time = 0;
float cycleTimeSec = 0;

// ************************************************************ REQUIRED INPUT FOR (pinv_damped FUNCTION) TO COMPUTE ANY DESIRED MATRIX INVERSE *************
double eps = 0.03;
int retPinvT;

//************************************************************ SETTING Kp AND Kd IN CARTESIAN SPACE ********************************************************
//float Kp = 300;
//float Kd = 15;
float Kpw = 700;
float Kdw = 100;

//************************************************************ PROOF OF ORTHOGHONAL FORCE VARIABLE *********************************************************
float dotProduct;

//************************************************************ SETTING Kp AND Kd IN JOINT SPACE ***********************************************************
float KpJ = 5;
float KdJ = 0.07;

//************************************************************ SETTING Kp AND Kd IN CARTESIAN SPACE WITH ORIENTATION CONTROL *******************************
float KpEU = 10;
float KdEU = 0.1;

//************************************************************ ROBOT REQUIRED VARIABLES SUCH AS ORCES, MATRICES AND ETC. **********************************
Matrix3d Mat;
MatrixXd trans_mat(6,6);
MatrixXd pinvT(3,3);
Vector3d eulerAngleOvershoot;
eulerAngleOvershoot << 0,0,0;
VectorXd svdT(3);
MatrixXd Jac(6,LBR_MNJ);
MatrixXd JacTr(LBR_MNJ,6);
MatrixXd JacEuler(6,LBR_MNJ);
MatrixXd JacEulerTr(LBR_MNJ,6);
MatrixXd JacEulerInv(LBR_MNJ,6);
MatrixXd T(3,3);
VectorXd sigma(6);
MatrixXd Rwb(3,3);
MatrixXd Fextb(6,1);
MatrixXd Fextw(6,1);
MatrixXd ExtTrq(7,1);
Vector3d PcWF;
Vector3d PrePcWF;
Vector3d P1WF;
Vector3d P2WF;
Vector3d P2primePreWF;
Vector3d P2primeWF;
Vector3d DeltaP;
Vector3d FpointDis;
MatrixXd ai(4,1);
MatrixXd ResMat(4,1);
MatrixXd TrajMat(4,4);
Vector3d PdesWF;
Vector3d PdesDotWF;
Vector3d Amp;
Vector3d Ave;
Vector3d P1TraWF;
Vector3d P2TraWF;
Vector3d qWF;
Vector3d qBF;
Vector3d Uw;
Vector3d Ub;
Vector3d Fw;
Vector3d Fb;
Vector3d Fdiff;
Vector3d FwLine;
Vector3d FbLine;
Vector3d Vec;
Vector3d VecEuler;
Vector3d PreVecEuler;
Vector3d CurVecEuler;
Vector3d DesVecEuler;
Vector3d PreVec;
Vector3d CurVec;
Vector3d DesVec;
Vector3d qCdot;
Vector3d qCcdot;
Vector3d qCcCdot;
Vector3d qCdotEuler;
MatrixXd DesTrq(7,1);
MatrixXd DesTrqEA(7,1);
MatrixXd F(6,1);
MatrixXd FEU(6,1);

float jointPosRad[LBR_MNJ];
float jointPosRadTminus1[LBR_MNJ];
double jointPosRad_d[LBR_MNJ];
double jointVelRads_d[LBR_MNJ];
double jointVelRads_dTminus1[LBR_MNJ];
double jointAccelRads2_d[LBR_MNJ];
float jointPosInitRad[LBR_MNJ];
float jointTrq[LBR_MNJ];
float jointExtTrqEst[LBR_MNJ];
float jointCmdPos[LBR_MNJ] = {0,0,0,0,0,0,0};
float jointCmdTrq[LBR_MNJ] = {0,0,0,0,0,0,0};
VectorXd jointValuesInit(LBR_MNJ);
//TIMER
timespec stampCurrentTime,stampInit,stampBegin,stampEnd,stampCycleDuration;
double cycleTimeNs = 0;
//PARAMETERS IDENTIFICATION
bool __initOk__ = false;
//JOINT SPEED LIMIT IN RAD/S
float jointVelocityLimit[LBR_MNJ];

jointVelocityLimit[0] = 1.92;
jointVelocityLimit[1] = 1.92;
jointVelocityLimit[2] = 2.23;
jointVelocityLimit[3] = 2.23;
jointVelocityLimit[4] = 3.56;
jointVelocityLimit[5] = 3.20;
jointVelocityLimit[6] = 3.20;

clock_gettime(CLOCK_REALTIME, &stampInit);

//************************************************************ ROBOT REQUIRED VARIABLES FOR LINE IMPLEMENTATION AND ADAPTATION *****************************

DeltaP(0) = 0;
DeltaP(1) = 0;
DeltaP(2) = 0;
P2WF(0) = 0.5;
P2WF(1) = 0.5;
P2WF(2) = 0.5;
P2primePreWF(0) = P2WF(0);
P2primePreWF(1) = P2WF(1);
P2primePreWF(2) = P2WF(2);

for(;;)
{
//GET THE STARTING TIME
clock_gettime(CLOCK_REALTIME, &stampCurrentTime);
stampBegin.tv_nsec = stampCurrentTime.tv_nsec - stampInit.tv_nsec;
stampBegin.tv_sec = stampCurrentTime.tv_sec - stampInit.tv_sec;

//************************************************************ PREVIOUS JOINT POSITIONS ********************************************************************
jointPrePosRadKuka[0] = jointPosRadKuka[0];
jointPrePosRadKuka[1] = jointPosRadKuka[1];
jointPrePosRadKuka[2] = jointPosRadKuka[2];
jointPrePosRadKuka[3] = jointPosRadKuka[3];
jointPrePosRadKuka[4] = jointPosRadKuka[4];
jointPrePosRadKuka[5] = jointPosRadKuka[5];
jointPrePosRadKuka[6] = jointPosRadKuka[6];

//************************************************************ ROTATION FROM WORLD FRAME TO BASE FRAME ******************************************************
Rwb(0,0) = 0.965925;
Rwb(0,1) = 0;
Rwb(0,2) = 0.258819;
Rwb(1,0) = 0;
Rwb(1,1) = 1;
Rwb(1,2) = 0;
Rwb(2,0) = -0.258819;
Rwb(2,1) = 0;
Rwb(2,2) = 0.965925;

//************************************************************ PREVIOUS END-EFFECTOR POSITION AND ORIENTATION ***********************************************
PrePcWF = Rwb * Vec;
PreVec = Rwb.inverse() * PrePcWF;
PreVecEuler = VecEuler;

//************************************************************ SYSTEM REAL TIME AND TIME CYCLE (sec) ********************************************************
Time = t;
cycleTimeSec = cycleTimeNs * 0.000000001;

kuka.doReceiveData();
kuka.displayComInfo();
if(kuka.getComQuality() == 3)
{
kuka.updateState();
kuka.getJointPos(jointPosRad);
kuka.getJointPos(jointPosRadKuka);
kuka.getJointTrq(jointTrq);
kuka.getJointExtTrqEst(jointExtTrqEst);

//************************************************************ EXTERNAL TORQUES READING FORM THE ROBOT *****************************************************
ExtTrq(0,0) = jointExtTrqEst[0];
ExtTrq(1,0) = jointExtTrqEst[1];
ExtTrq(2,0) = jointExtTrqEst[2];
ExtTrq(3,0) = -jointExtTrqEst[3];
ExtTrq(4,0) = jointExtTrqEst[4];
ExtTrq(5,0) = jointExtTrqEst[5];
ExtTrq(6,0) = jointExtTrqEst[6];

//************************************************************ CALCULATION OF JACOBIAN TRANSPOSE AND JACOBIAN EULER TRANSPOSE *******************************
kuka.directKinematics(Mat,Vec);
kuka.computeEulerAngles(Mat,VecEuler);
kuka.jacobian(Jac);
JacTr = Jac.transpose();
T << cos(VecEuler(1))*cos(VecEuler(2)), -sin(VecEuler(2)), 0,
cos(VecEuler(1))*sin(VecEuler(2)), cos(VecEuler(2)), 0,
-sin(VecEuler(1)),
0,
1;
retPinvT = pinv_damped(T, pinvT, svdT, eps);
trans_mat.topLeftCorner(3,3) = MatrixXd::Identity(3,3);
trans_mat.topRightCorner(3,3) = MatrixXd::Zero(3,3);
trans_mat.bottomLeftCorner(3,3) = MatrixXd::Zero(3,3);
trans_mat.bottomRightCorner(3,3) = pinvT;
JacEuler = trans_mat * Jac;
JacEulerTr = JacEuler.transpose();
pinv_damped(JacEulerTr, JacEulerInv, sigma, eps);

//************************************************************ ECALCULATION OF EXTERNAL FORCES **************************************************************
Fextb = JacEulerInv * ExtTrq;
Fextb(0,0) = - Fextb(0,0);
Fextb(1,0) = - Fextb(1,0);
Fextb(2,0) = - Fextb(2,0);

//************************************************************ LINE ADAPTATION USING HUMAN-INTERACTION-FORCE ************************************************
PcWF = Rwb * Vec;
P1WF(0) = -0.5;
P1WF(1) = 0.5;
P1WF(2) = 0.5;
float FT_limit = 5;
float FT_gain = 0.0001;
if(fabs(Fextb(0,0)) > FT_limit){
DeltaP(0) = FT_gain * Fextb(0,0);
}
else
DeltaP(0) = 0;
if(fabs(Fextb(1,0)) > FT_limit){
    DeltaP(1) = FT_gain * Fextb(1,0);
}
else
DeltaP(1) = 0;
if(fabs(Fextb(2,0)) > FT_limit){
DeltaP(2) = FT_gain * Fextb(2,0);
}
else
DeltaP(2) = 0;
/// ######################################## ADAPTATION PROCESS ##################################################
P2primeWF(0) = P2primePreWF(0) + DeltaP(0);
P2primeWF(1) = P2primePreWF(1) + DeltaP(1);
P2primeWF(2) = P2primePreWF(2) + DeltaP(2);
/// ###########################################################################################################
if(P2primeWF(0) > 0.5){
P2primeWF(0) = 0.5;
}
if(P2primeWF(0) < 0){
P2primeWF(0) = 0;
}
if(P2primeWF(1) > 0.6){
P2primeWF(1) = 0.6;
}
if(P2primeWF(1) < 0.3){
P2primeWF(1) = 0.3;
}
if(P2primeWF(2) > 0.7){
P2primeWF(2) = 0.7;
}
if(P2primeWF(2) < 0.3){
P2primeWF(2) = 0.3;
}

/************************************************************ EXAMPLES OF MOVING ROBOT ALONG
AXES **********************************************************
//******* MOVING ALONG Y-AXIS
/*
P1WF(0) = 0;
P1WF(1) = -0.5;
P1WF(2) = 0.5;
P2WF(0) = 0;
P2WF(1) = 0.5;
P2WF(2) = 0.5;
*/
//****************************
//******* MOVING ALONG Z-AXIS
/*
P1WF(0) = 0;
P1WF(1) = 0.5;
P1WF(2) = -0.5;
P2WF(0) = 0;
P2WF(1) = 0.5;
P2WF(2) = 0.5;
*/
//****************************

//******* MOVING ALONG 2D-LINE (YZ PLANE)
/*
P1WF(0) = 0;
P1WF(1) = 0;
P1WF(2) = 0;
P2WF(0) = 0;
P2WF(1) = 0.5;
P2WF(2) = 0.5;
*/
//****************************
//******* MOVING ALONG ARBITRAY LINE IN 3D
/*
P1WF(0) = -0.5;
P1WF(1) = -0.5;
P1WF(2) = -0.5;
P2WF(0) = 0.5;
P2WF(1) = 0.5;
P2WF(2) = 0.5;
*/
//****************************

//************************************************************ PROJECTION POINT CALCULATION IN WORLD FRAME **************************************************
qWF(0) = (P1WF(0)*pow(P2primeWF(2),2) + ((-P2primeWF(0) - P1WF(0))*P1WF(2) + (P2primeWF(0)
- P1WF(0))*PcWF(2))*P2primeWF(2) + P2primeWF(0)*pow(P1WF(2),2) + (P1WF(0) - P2primeWF(0))*PcWF(2)*P1WF
(2) + P1WF(0)*pow(P2primeWF(1),2) + ((-P2primeWF(0) - P1WF(0))*P1WF(1) + (P2primeWF(0) - P1WF(0))*PcWF
(1))*P2primeWF(1) + P2primeWF(0)*pow(P1WF(1),2) + (P1WF(0) - P2primeWF(0))*PcWF(1)*P1WF(1) + PcWF
(0)*pow(P2primeWF(0),2) - 2*PcWF(0)*P1WF(0)*P2primeWF(0) + PcWF(0)*pow(P1WF(0),2))/
(pow(P2primeWF(2),2) - 2*P1WF(2)*P2primeWF(2) + pow(P1WF(2),2)+pow(P2primeWF(1),2) - 2*P1WF
(1)*P2primeWF(1) + pow(P1WF(1),2) + pow(P2primeWF(0),2) - 2*P1WF(0)*P2primeWF(0) + pow(P1WF(0),2));

qWF(1) = (P1WF(1)*pow(P2primeWF(2),2) + ((-P2primeWF(1) - P1WF(1))*P1WF(2) + (P2primeWF(1)
- P1WF(1))*PcWF(2))*P2primeWF(2) + P2primeWF(1)*pow(P1WF(2),2) + (P1WF(1) - P2primeWF(1))*PcWF(2)*P1WF
(2) +PcWF(1)*pow(P2primeWF(1),2) + (-2*PcWF(1)*P1WF(1) + (PcWF(0) - P1WF(0)*P2primeWF(0)) + pow(P1WF
(0),2) - PcWF(0)*P1WF(0))*P2primeWF(1) + PcWF(1)*pow(P1WF(1),2) + (pow(P2primeWF(0),2) + (-P1WF(0) -
PcWF(0))*P2primeWF(0) + PcWF(0)*P1WF(0))*P1WF(1))/
(pow(P2primeWF(2),2) - 2*P1WF(2)*P2primeWF(2) + pow(P1WF(2),2)+pow(P2primeWF(1),2) - 2*P1WF
(1)*P2primeWF(1) + pow(P1WF(1),2) + pow(P2primeWF(0),2) - 2*P1WF(0)*P2primeWF(0) + pow(P1WF(0),2));

qWF(2) = (PcWF(2)*pow(P2primeWF(2),2) + (-2*PcWF(2)*P1WF(2) + (PcWF(1) - P1WF(1))*P2primeWF
(1) + pow(P1WF(1),2) - PcWF(1)*P1WF(1) + (PcWF(0) - P1WF(0))*P2primeWF(0) + pow(P1WF(0),2) - PcWF
(0)*P1WF(0))*P2primeWF(2) + PcWF(2)*pow(P1WF(2),2) + (pow(P2primeWF(1),2) + (-P1WF(1) - PcWF(1)) *
P2primeWF(1) + PcWF(1)*P1WF(1) + pow(P2primeWF(0),2) + (-P1WF(0) - PcWF(0))*P2primeWF(0) + PcWF(0)*P1WF
(0))*P1WF(2))/
(pow(P2primeWF(2),2) - 2*P1WF(2)*P2primeWF(2) + pow(P1WF(2),2)+pow(P2primeWF(1),2) - 2*P1WF
(1)*P2primeWF(1) + pow(P1WF(1),2) + pow(P2primeWF(0),2) - 2*P1WF(0)*P2primeWF(0) + pow(P1WF(0),2));

//************************************************************ PROJECTION POINT CALCULATION IN BASE FRAME ***************************************************
qBF = Rwb.inverse() * qWF;

//************************************************************ REQUIRED FORCE FOR REACHING THE LINE AND STAYING ON IT ***************************************
Fw = Kpw * (qWF - PcWF) - Kdw * ((PcWF - PrePcWF)/cycleTimeSec);

//************************************************************ PROOF OF ORTHOGONAL FORCE TO THE LINE ********************************************************
dotProduct = Fw(0) * (P2WF(0) - P1WF(0)) + Fw(1) * (P2WF(1) - P1WF(1)) + Fw(2) * (P2WF(2) - P1WF(2));

//************************************************************ IMPLEMENTING LINEAR TRAJECTORY ALONG X-AXIS **************************************************
/*
if (Time >= 0 and Time <= 5){
ai(0)
ai(1)
ai(2)
ai(3)
=
=
=
=
-0.5;
0;
0.06;
-0.008;
TrajMat(0,0)
TrajMat(0,1)
TrajMat(0,2)
TrajMat(0,3) =
=
=
= 1;
0;
0;
0;
TrajMat(1,0)
TrajMat(1,1)
TrajMat(1,2)
TrajMat(1,3) =
=
=
= 0;
1;
0;
0;
TrajMat(2,0)
TrajMat(2,1)
TrajMat(2,2)
TrajMat(2,3) =
=
=
= 1;
Time;
pow(Time,2);
pow(Time,3);
TrajMat(3,0)
TrajMat(3,1)
TrajMat(3,2)
TrajMat(3,3) =
=
=
= 0;
1;
2 * Time;
3 * pow(Time,2);
ResMat = TrajMat * ai;
FwLine(0) =
//FwLine(0)
FwLine(1) =
FwLine(2) =
}
50 * (ResMat(2) - qWF(0));
= 25 * (ResMat(3));
0;
0;
else {
FwLine(0) = 0;
FwLine(1) = 0;
FwLine(2) = 0;
}
*/

//************************************************************ IMPLEMENTING LINEAR PERIODIC TRAJECTORY IN 3D ************************************************
Amp(0) = P2primeWF(0) - P1WF(0);
Amp(1) = P2primeWF(1) - P1WF(1);
Amp(2) = P2primeWF(2) - P1WF(2);
Ave(0) = (P2primeWF(0) + P1WF(0))/2;
Ave(1) = (P2primeWF(1) + P1WF(1))/2;
Ave(2) = (P2primeWF(2) + P1WF(2))/2;
PdesWF(0) = Amp(0) * 0.5 * sin(0.35 * Time) + Ave(0);
PdesDotWF(0) = 0.35 * 0.5 * Amp(0) * cos(0.75 * Time);
PdesWF(1) = Amp(1) * 0.5 * sin(0.35 * Time) + Ave(1);
PdesDotWF(1) = 0.35 * 0.5 * Amp(1) * cos(0.75 * Time);
PdesWF(2) = Amp(2) * 0.5 * sin(0.35 * Time) + Ave(2);
PdesDotWF(2) = 0.35 * 0.5 * Amp(2) * cos(0.75 * Time);
// ******* NOISE REMOVAL
if(PdesWF(0) != 0 && fabs(PdesWF(1) - PcWF(1) > 0.5) && fabs(PdesWF(2) - PcWF(2) > 0.5)){
PdesWF(0) = PdesWF(0);
PdesWF(1) = 0;
PdesWF(2) = 0;
}
if(fabs(PdesWF(0) - PcWF(0) > 0.5) && PdesWF(1) != 0 && fabs(PdesWF(2) - PcWF(2) > 0.5)){
PdesWF(0) = 0;
PdesWF(1) = PdesWF(1);
PdesWF(2) = 0;
}
if(fabs(PdesWF(0) - PcWF(0) > 0.5) && fabs(PdesWF(1) - PcWF(1) > 0.5) && PdesWF(2) != 0){
PdesWF(0) = 0;
PdesWF(1) = 0;
PdesWF(2) = PdesWF(2);
}
if(PdesWF(0) !=
PdesWF(0) =
PdesWF(1) =
PdesWF(2) =
} 0 && PdesWF(1) != 0 && fabs(PdesWF(2) - PcWF(2) > 0.5)){
PdesWF(0);
PdesWF(1);
0;
if(PdesWF(0) !=
PdesWF(0) =
PdesWF(1) =
PdesWF(2) =
} 0 && fabs(PdesWF(1) - PcWF(1) > 0.5) && PdesWF(2) != 0){
PdesWF(0);
0;
PdesWF(2);
if(fabs(PdesWF(0) - PcWF(0) > 0.5) && PdesWF(1) != 0 && PdesWF(2) != 0){
PdesWF(0) = 0;
PdesWF(1) = PdesWF(1);
PdesWF(2) = PdesWF(2);
}
if(PdesWF(0) != 0 && PdesWF(1) != 0 && PdesWF(2) != 0){
PdesWF(0) = PdesWF(0);
PdesWF(1) = PdesWF(1);
PdesWF(2) = PdesWF(2);
}

// ******* REQUIRED FORCE FOR IMPLEMENTING THE PERIODIC TRAJECTORY ALONG AN ARBITRAY LINE IN 3D IN WORLD FRAME
FwLine(0) = 700 * (PdesWF(0) - qWF(0)) + 100 * (PdesDotWF(0) - ((PcWF(0) - PrePcWF(0))/
cycleTimeSec));
FwLine(1) = 700 * (PdesWF(1) - qWF(1)) + 100 * (PdesDotWF(1) - ((PcWF(1) - PrePcWF(1))/
cycleTimeSec));
FwLine(2) = 700 * (PdesWF(2) - qWF(2)) + 100 * (PdesDotWF(2) - ((PcWF(2) - PrePcWF(2))/
cycleTimeSec));
//
************************************************************************************************
if(PdesWF (0) == 0){
FwLine(0) = 0;
}
if(PdesWF (1) == 0){
FwLine(1) = 0;
}
if(PdesWF (2) == 0){
FwLine(2) = 0;
}
// ******* REQUIRED FORCE FOR IMPLEMENTING THE PERIODIC TRAJECTORY ALONG AN ARBITRAY LINE
IN 3D IN BASE FRAME
FbLine = Rwb.inverse() * FwLine;
//
**********************************************************************************************************
// ******* DIFFERENCE BETWEEN EXTERNAL FORCE AND LINE REQUIRED FORCE
Fdiff(0) = Fextb(0,0) - Fb(0);
Fdiff(1) = Fextb(1,0) - Fb(1);
Fdiff(2) = Fextb(2,0) - Fb(2);
// *****************************************************************
// ******* REQUIRED LINE FORCE CALCULATION IN BASE FRAME
Fb = Rwb.inverse() * Fw;

///************************************************************ MY MINI-PROJECT PART (NOT
REQUIRED HERE - ORIENTATION CONTROL AND APPLYING GRAVITY ALONG ANY DESIRED LINE) ***
//******* ORIENTATION CONTROL (ROBOT BASE FRAME)*******
//CurVecEuler(0) = VecEuler(0);
//CurVecEuler(1) = VecEuler(1);
//CurVecEuler(2) = VecEuler(2);
//DesVecEuler(0) = -1.767780;
//DesVecEuler(1) = 1.237793;
//DesVecEuler(2) = -0.211584;
//************ CURRENT END-EFFECTOR VELOCITY (qCdot CARTESIAN SPACE) ************
//qCdot(0) = (CurVec(0) - PreVec(0)) / cycleTimeSec;
//qCdot(1) = (CurVec(1) - PreVec(1)) / cycleTimeSec;
//qCdot(2) = (CurVec(2) - PreVec(2)) / cycleTimeSec;
/*
if(Time >= 0){
qCdot(0) = -0.292914;
//(CurVec(0) - PreVec(0)) / cycleTimeSec;
qCdot(1) = 0;
//(CurVec(1) - PreVec(1)) / cycleTimeSec;
qCdot(2) = 0;
//(CurVec(2) - PreVec(2)) / cycleTimeSec;
}

if(Time >= 10){
qCdot(0) = 0;
//(CurVec(0) - PreVec(0)) / cycleTimeSec;
qCdot(1) = 0;
//(CurVec(1) - PreVec(1)) / cycleTimeSec;
qCdot(2) = 0.297082;
//(CurVec(2) - PreVec(2)) / cycleTimeSec;
}

if(Time >= 20){
qCdot(0) = -0.282914;
//(CurVec(0) - PreVec(0)) / cycleTimeSec;
qCdot(1) = 0;
//(CurVec(1) - PreVec(1)) / cycleTimeSec;
qCdot(2) = 0;
//(CurVec(2) - PreVec(2)) / cycleTimeSec;
}

if(Time >= 30){cartesian interpolated trajectory
qCdot(0) = 0;
//(CurVec(0) - PreVec(0)) / cycleTimeSec;
qCdot(1) = 0;
//(CurVec(1) - PreVec(1)) / cycleTimeSec;
qCdot(2) = 0.287082;
//(CurVec(2) - PreVec(2)) / cycleTimeSec;
}

if(Time >= 40){
qCdot(0) = (CurVec(0) - PreVec(0)) / cycleTimeSec;
qCdot(1) = (CurVec(1) - PreVec(1)) / cycleTimeSec;
qCdot(2) = (CurVec(2cartesian interpolated trajectory) - PreVec(2)) / cycleTimeSec;
}

*/
//qCcdot(0) = (CurVecEuler(0) - PreVecEuler(0)) / cycleTimeSec;
//qCcdot(1) = (CurVecEuler(1) - PreVecEuler(1)) / cycleTimeSec;
//qCcdot(2) = (CurVecEuler(2) - PreVecEuler(2)) / cycleTimeSec;
//qCcCdot(0) = (DesVecEuler(0) - CurVecEuler(0)) / cycleTimeSec;
//qCcCdot(1) = (DesVecEuler(1) - CurVecEuler(1)) / cycleTimeSec;
//qCcCdot(2) = (DesVecEuler(2) - CurVecEuler(2)) / cycleTimeSec;
//***************************************************************

//********** MOVING AROUND THE DESIRED POINT *******
//F(0,0) = Kp * (DesVec(0) - CurVec(0)) + Kd * (0 - qCdot(0));
//F(1,0) = Kp * (DesVec(1) - CurVec(1)) + Kd * (0 - qCdot(1));
//F(2,0) = Kp * (DesVec(2) - CurVec(2)) + Kd * (0 - qCdot(2));

//********** MOVING ALONG X-AXIS *******
//F(0,0) = 0;
//Kp * (DesVec(0) - CurVec(0)) + Kd * (0 - qCdot(0));
//F(1,0) = Kp * (DesVec(1) - CurVec(1)) + Kd * (0 - qCdot(1));
//F(2,0) = Kp * (DesVec(2) - CurVec(2)) + Kd * (0 - qCdot(2));
//***************************************************************
//********** MOVING ALONG Y-AXIS *******
//F(0,0) = Kp * (DesVec(0) - CurVec(0)) + Kd * (0 - qCdot(0));
//F(1,0) = 0;
//Kp * (DesVec(1) - CurVec(1)) + Kd * (0 - qCdot(1));
//F(2,0) = Kp * (DesVec(2) - CurVec(2)) + Kd * (0 - qCdot(2));
//***************************************************************
//********** MOVING ALONG Z-AXIS *******
//F(0,0) = Kp * (DesVec(0) - CurVec(0)) + Kd * (0 - qCdot(0));
//F(1,0) = Kp * (DesVec(1) - CurVec(1)) + Kd * (0 - qCdot(1));
//F(2,0) = 0;
//Kp * (DesVec(2) - CurVec(2)) + Kd * (0 - qCdot(2));
//***************************************************************
//********** SETTING THE END-EFFECTOR ORIENTATION NOT FIXED **********
//F(3,0) = 0;
//F(4,0) = 0;
//F(5,0) = 0;
//********** THE END-EFFECTOR POSITION AND ORIENTATION CONTROL **********
//FEU(0,0) = Kp * (DesVec(0) - CurVec(0)) + Kd * (qCcCdot(0) - qCdot(0));
//FEU(1,0) = Kp * (DesVec(1) - CurVec(1)) + Kd * (qCcCdot(1) - qCdot(1));
//FEU(2,0) = Kp * (DesVec(2) - CurVec(2)) + Kd * (qCcCdot(2) - qCdot(2));

//************************************************************ FORCE CONSTRAINTS AND FORCE CALCULATION *****************************************************
if(fabs(Fextb(0,0)) > FT_limit){
FbLine(0) = 0;
}
if(fabs(Fextb(1,0)) > FT_limit){
FbLine(1) = 0;
}
if(fabs(Fextb(2,0)) > FT_limit){
FbLine(2) = 0;
}

///######################################### PERIODIC TRAJECTORY ##################################################
/*
FbLine(0) = 0;
FbLine(1) = 0;
FbLine(2) = 0;
*/
///###########################################################################################################
FEU(0,0) = Fb(0) + FbLine(0);
FEU(1,0) = Fb(1) + FbLine(1);
FEU(2,0) = Fb(2) + FbLine(2);
//***********************************************************************************************************
///************************************************************ MY MINI-PROJECT PART (NOT REQUIRED HERE - ORIENTATION CONTROL AND APPLYING GRAVITY ALONG ANY DESIRED LINE) ***
//FEU(0,0) = Kp * (DesVecBF(0) - CurVecBF(0)) + Kd * (0 - qCdot(0));
//FEU(1,0) = Kp * (DesVecBF(1) - CurVecBF(1)) + Kd * (0 - qCdot(1));
//FEU(2,0) = Kp * (DesVecBF(2) - CurVecBF(2)) + Kd * (0 - qCdot(2));
//FEU(0,0) = Kp * (DesVec(0) - CurVec(0)) + Kd * (0 - qCdot(0));
//FEU(1,0) = Kp * (DesVec(1) - CurVec(1)) + Kd * (0 - qCdot(1));
//FEU(2,0) = Kp * (DesVec(2) - CurVec(2)) + Kd * (0 - qCdot(2));

FEU(3,0) = 0;//KpEU * (DesVecEuler(0) - CurVecEuler(0)) + KdEU * (0 - qCcdot(0));
FEU(4,0) = 0;//KpEU * (DesVecEuler(1) - CurVecEuler(1)) + KdEU * (0 - qCcdot(1));
FEU(5,0) = 0;//KpEU * (DesVecEuler(2) - CurVecEuler(2)) + KdEU * (0 - qCcdot(2));

float lim[3] = {M_PI,100000,M_PI};
for(int i = 0; i < 3; i++)
{
if(fabs((CurVecEuler(i) - PreVecEuler(i))) > lim[i])
{
if((CurVecEuler(i) - PreVecEuler(i)) > 0)
{
//cout << endl << "OVERSHOOT! + axis" << i <<endl;
eulerAngleOvershoot(i) = eulerAngleOvershoot(i) + 1;
}
else
{
//cout << endl << "OVERSHOOT! - axis" << i << endl;
eulerAngleOvershoot(i) = eulerAngleOvershoot(i) - 1;
}
}
}
jointPosRadKuka[3] = -jointPosRadKuka[3];
if(!__initOk__)
{
for(int i=0;i<LBR_MNJ;i++)
{
jointPosInitRad[i] = jointPosRad[i];
jointValuesInit(i) = jointPosRad[i];
}
if(kuka_mode == POSITION)
{
//POSITION CONTROL
for(int i=0;i<LBR_MNJ;i++)
jointCmdPos[i] = jointPosRad[i];
jointCmdPos[3] = -jointCmdPos[3];
kuka.doPositionControl(jointCmdPos);
if(t > 5)
{
__initOk__ = true;
t = 0;
}
}
else //TORQUE CONTROL
{
if(kuka.getCurrentControlStrategy() != JOINT_IMPEDANCE_CONTROL)
{
kuka.doChangeControlStrategy(JOINT_IMPEDANCE_CONTROL);
t = 0;
}
else
{
kuka.doJntImpedanceControl(jointPosRadKuka, NULL, NULL, jointCmdTrq);
if(t > 2)
{
__initOk__ = true;
t = 0;
}
}
}
}
else
{
    //************************************************************ CALCULATION OF DESIRED JOINT TORQUES TO COMPLETE THE TASK ***********************************
DesTrqEA = JacEulerTr * FEU;
jointDesTrq[0] = DesTrqEA(0);
jointDesTrq[1] = DesTrqEA(1);
jointDesTrq[2] = DesTrqEA(2);
jointDesTrq[3] = DesTrqEA(3);
jointDesTrq[4] = DesTrqEA(4);
jointDesTrq[5] = DesTrqEA(5);
jointDesTrq[6] = DesTrqEA(6);

for(int i = 0; i < 7 ;i++)
{
if(fabs(jointDesTrq[i]) > TORQUELIMIT)
{
if(jointDesTrq[i] > 0)
jointCmdTrq[i] = TORQUELIMIT;
else
jointCmdTrq[i] = -TORQUELIMIT;
}
else
{
jointCmdTrq[i] = jointDesTrq[i];
}
jointCmdPos[i] = jointPosRad[i] + 0.0;
}
/** Send the torque command*/
//REMEMBER THAT THETA4 IS INVERTED
jointCmdTrq[3] = -jointCmdTrq[3] ;
jointCmdPos[3] = -jointCmdPos[3];
if(kuka_mode == POSITION)
{
//CONTROL THE KUKA ROBOT
if(kuka.getCurrentControlStrategy() != POSITION_CONTROL)
{
kuka.doChangeControlStrategy(POSITION_CONTROL);
}
else
{
kuka.doPositionControl(jointCmdPos);
}
}
else //POSITION CONTROL
{
//CONTROL THE KUKA ROBOT
if(kuka.getCurrentControlStrategy() != JOINT_IMPEDANCE_CONTROL)
{
kuka.doChangeControlStrategy(JOINT_IMPEDANCE_CONTROL);
}
else
{
kuka.doJntImpedanceControl(jointPosRadKuka, NULL, NULL, jointCmdTrq);
}
}
}
//************************************************************ CALCULATION OF PREVIOUSPOINT ***************************************************************
P2primePreWF(0) = P2primeWF(0);
P2primePreWF(1) = P2primeWF(1);
P2primePreWF(2) = P2primeWF(2);
}

printf(" Time %lf P2primeX %lf P2primeY %lf P2primeZ %lf\n",Time,P2primeWF(0),P2primeWF
(1),P2primeWF(2));
printf(" Time %lf PdesX %lf PdesY %lf PdesZ %lf PcWFX %lf PcWFY %lf PcWFZ %lf\n",Time,PdesWF
(0),PdesWF(1),PdesWF(2),PcWF(0),PcWF(1),PcWF(2));
printf(" Time %lf FbX %lf FbY %lf FbZ %lf FextX %lf FextY %lf FextZ %lf\n",Time,Fb(0),Fb(1),Fb
(2),Fextb(0,0),Fextb(1,0),Fextb(2,0));

//SEND DATA TO KRL
kuka.doSendData();
//END TIME
clock_gettime(CLOCK_REALTIME, &stampCurrentTime);
stampEnd.tv_nsec = stampCurrentTime.tv_nsec - stampInit.tv_nsec;
stampEnd.tv_sec = stampCurrentTime.tv_sec - stampInit.tv_sec;
//CYCLE DURATION
if ((stampEnd.tv_sec - stampBegin.tv_sec) == 0)
stampCycleDuration.tv_nsec = stampEnd.tv_nsec - stampBegin.tv_nsec;
else
stampCycleDuration.tv_nsec = stampEnd.tv_nsec - stampBegin.tv_nsec + 1000000000;
cycleTimeNs = stampCycleDuration.tv_nsec;
tminus1 = t;
t = t + cycleTimeNs*0.000000001;
}
return EXIT_SUCCESS;
}

//********** IT WILL BE USED INSTEAD OF THE INVERSE FUNCTION **********
bool pinv_damped(const MatrixXd A, MatrixXd &invA, VectorXd &sigma, double eps)
{
MatrixXd toInvert;
double lambda;
int ret;
int m;
Eigen::JacobiSVD<Eigen::MatrixXd> svd_A;
svd_A.compute(A.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
sigma = svd_A.singularValues();
//sigma(m);
if(A.rows() <= A.cols())
{
m = (int)A.rows() - 1;
toInvert = A*A.transpose();
}
else
{
m = (int)A.cols() - 1;
toInvert = A.transpose()*A;
}
if (sigma(m) <= eps)
{
lambda = pow(eps,2) - pow(sigma(m),2);
ret = true;
}
else
{
lambda = 0;
ret = false;
}
MatrixXd dampingMat = pow(lambda,2)*MatrixXd::Identity(6,6);
for(int i=0;i< toInvert.rows();i++)
toInvert(i,i) = toInvert(i,i) + pow(lambda,2);
if(A.rows() <= A.cols())
invA = A.transpose()*(toInvert.inverse());
else
invA = (toInvert.inverse())*A.transpose();
return ret;
}