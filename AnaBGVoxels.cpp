/*********************************************************************
 * Author           : Lan-sx
 * Email            : shexin@ihep.ac.cn
 * Created          : 2025-03-21 23:20
 * Filename         : AnaBGVoxels.cpp
 * Description      : 
 * ******************************************************************/
//std
#include <vector>
#include <array>
#include <tuple>
#include <map>
#include <memory>
#include <algorithm>

//ROOT CERN
#include "TApplication.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGeoManager.h"
#include "TEnv.h"
#include "TPCVoxelHelper.hh"
#include "TGLTH3Composition.h"
#include "TDatime.h"
#include "TGLAxis.h"
#include "Lansxlogon.h"

#include "PlotVoxelsGL.hh"
using namespace std;

int main(int argc, char **argv)
{
    LansxFormat::myStyle();
    TColor::InvertPalette();
    gEnv->SetValue("OpenGL.SelectionBufferSize", 4096);
    auto file = TFile::Open("../../Higgs/MixedTestPolar_1mm6mm2mm_NearSearch026_Trandom.root");
    if(!file || file->IsZombie())
    {
        std::printf("[error] file does not exist!\n");
        return -1;
        //return;
    }
    auto tsig = dynamic_cast<TTree*>(file->Get("sigbgVoxels"));
    auto tbkg = dynamic_cast<TNtuple*>(file->Get("tbkg"));
    if(!tsig || !tbkg)
    {
        std::printf("[error] Tree sigbgVoxels does not exist\n");
        return -1;
        //return ;
    }

    TApplication app("app",&argc, argv);
    auto helper = TPCVoxelHelper::GetInstance();
    helper->SetTPCVoxel(1, 0.6, 0.1);
    std::printf("[info] pad type:%d, padlength=%.4f, padwidth=%.4f\n", helper->GetVoxelType(),
        helper->GetPadLenth(),
        helper->GetPadWidth());
   
    gStyle->SetCanvasPreferGL(kTRUE);
   
    auto myc = new TCanvas("myc", "myc", 800, 800);
    myc->SetFillColor(kBlack);

    PlotVoxelsGL pltVoxelgl(tsig, tbkg);
    pltVoxelgl.EnableSigVoxelColorPalette();
    //pltVoxelgl.AddAnnotation("test annotation");
    pltVoxelgl.SetPhiRange(5.);
    pltVoxelgl.Enable3D();
    auto topnode = pltVoxelgl.ViewVoxels(1000);
    topnode->Draw("ogl");
    auto glview = dynamic_cast<TGLViewer*>(myc->GetViewer3D());
    pltVoxelgl.SetGLviewer(glview);
    TDatime time;
    pltVoxelgl.AddAnnotation(time.AsString());
    //glview->SavePicture("test.png");
    //glview->SetStyle(TGLRnrCtx::kWireFrame);
    //glview->SetStyle(TGLRnrCtx::kOutline);

    std::printf("[info] ====> Code finished Press Ctrl+C to exit :)");
    app.Run(1);
    return 0;
}
