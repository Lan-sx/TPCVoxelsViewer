/*********************************************************************
 * Author           : Lan-sx
 * Email            : shexin@ihep.ac.cn
 * Created          : 2025-04-12 00:42
 * Filename         : PlotVoxelsGL.hh
 * Description      : 
 * Update           : 
 * ******************************************************************/
#ifndef __PLOTVOXELSGL__
#define  __PLOTVOXELSGL__ 1
//std
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

//ROOT CERN
#include "TStyle.h"
#include "TMath.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TNode.h"
#include "TTUBS.h"
#include "TBRIK.h"
#include "TGLViewer.h"
#include "TGLAnnotation.h"

//Users
#include "TPCVoxelHelper.hh"

class PlotVoxelsGL
{
public:
    PlotVoxelsGL(TTree* tsig, TNtuple* tbkg);
    ~PlotVoxelsGL() {}

    //public methods
    void SetGLviewer(TGLViewer* glviewer) { m_glviewer = glviewer; }
    void Enable3D() { m_enable3D = true; }
    void EnableSigVoxelColorPalette() { m_sigVoxelPal = true; }
    void SetPhiRange(double phi) { m_phiRange = phi; }
    void AddAnnotation(const char* txt);
    TNode* ViewVoxels(int trkid=0);

private:
    TGLViewer* m_glviewer;
    TTree* m_tsig;
    TNtuple* m_tbkg;
    bool m_enable3D, m_sigVoxelPal;
    double m_phiRange;

};

#endif
