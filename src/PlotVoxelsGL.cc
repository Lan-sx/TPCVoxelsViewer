/*********************************************************************
 * Author           : Lan-sx
 * Email            : shexin@ihep.ac.cn
 * Created          : 2025-04-12 00:42
 * Filename         : PlotVoxelsGL.cpp
 * Description      : 
 * Update           : 
 * ******************************************************************/
#include "PlotVoxelsGL.hh"

PlotVoxelsGL::PlotVoxelsGL(TTree* tsig, TNtuple* tbkg) : m_glviewer(nullptr), m_tsig(tsig), m_tbkg(tbkg)
{
    m_enable3D = false;
    m_sigVoxelPal = false;
    m_phiRange = 5.;
}

TNode* PlotVoxelsGL::ViewVoxels(int trk)
{
    int trkIdx(0);
    double trkMomentum(0.), trkTheta(0.), trkPhi(0.);
    auto vPixelIdx = new std::vector<int>();
    auto vZbinIdx = new std::vector<int>();
    auto vSigQ = new std::vector<int>();

    m_tsig->SetBranchAddress("trkIdx", &trkIdx);
    m_tsig->SetBranchAddress("trkMomentum", &trkMomentum);
    m_tsig->SetBranchAddress("trkTheta", &trkTheta);
    m_tsig->SetBranchAddress("trkPhi", &trkPhi);
    m_tsig->SetBranchAddress("PixelIdx", &vPixelIdx);
    m_tsig->SetBranchAddress("ZbinIdx", &vZbinIdx);
    m_tsig->SetBranchAddress("SigQ", &vSigQ);

    float m_pixelIdx(0.), m_zbinIdx(0.), m_charge(0.);
    m_tbkg->SetBranchAddress("pixelIdx/F", &m_pixelIdx);
    m_tbkg->SetBranchAddress("zbinIdx/F", &m_zbinIdx);
    m_tbkg->SetBranchAddress("charge/F", &m_charge);

    auto helper = TPCVoxelHelper::GetInstance();
    auto padtype = helper->GetVoxelType();

    m_tsig->GetEntry(trk);
    auto MaxQ = std::max_element(vSigQ->begin(), vSigQ->end());
    //auto xyidx0 = vPixelIdx->at(0);
    //auto zidx0 = vZbinIdx->at(0);

    //auto padvertex0 = helper->GetPadVertex(xyidx0);
    std::printf("[info] Trk Theta=%.4f, Phi=%.4f\n", trkTheta, trkPhi);
    //std::printf("[info] r1=%.2f,r2=%.2f,phi1=%.4f,phi2=%.4f\n", std::get<0>(padvertex0),
    //    std::get<1>(padvertex0),
    //    std::get<2>(padvertex0) * TMath::RadToDeg(),
    //    std::get<3>(padvertex0) * TMath::RadToDeg());

    auto topShape = new TTUBS("topShape", "void", "vacuum", _TPCinnerR, _TPCouterR, 500, 0., 360.);
    auto topNode = new TNode("topNode", "", "topShape", 0, 0, 0);
    topNode->cd();
    std::vector<int> vzbinSig;
    
    double x_node(0.), y_node(0.), z_node(0.);

    //sig
    for (size_t ii = 0; ii < vPixelIdx->size(); ++ii)
    {
        int xyidx0 = vPixelIdx->at(ii);
        int zidx0 = vZbinIdx->at(ii);
        auto charge_ii = vSigQ->at(ii);
        if (xyidx0 >= 0 && zidx0 >= 0)
        {
            auto padvertex_ii = helper->GetPadVertex(xyidx0);
            double r_start = std::get<0>(padvertex_ii);
            double r_end = std::get<1>(padvertex_ii);
            double phi_start = std::get<2>(padvertex_ii) * TMath::RadToDeg();
            double phi_end = std::get<3>(padvertex_ii) * TMath::RadToDeg();
            auto linecolor = m_sigVoxelPal ? gStyle->GetColorPalette( 256 *  (double(charge_ii) / (*MaxQ)) ) : kMagenta;
            //auto linecolor = m_sigVoxelPal ? gStyle->GetColorPalette(nColors*100) : kMagenta;
            if (padtype == _Polar)
            {
                auto voxels_ii = new TTUBS(Form("voxels_%d", ii), "void", "vacuum", r_start, r_end, 0.5, phi_start, phi_end);
                voxels_ii->SetLineColor(linecolor);
            }
            else if (padtype == _Cartesian)
            {
                auto voxels_ii = new TBRIK(Form("voxels_%d", ii), "void", "vacuum", helper->GetPadWidth() * 0.5, helper->GetPadWidth() * 0.5, 0.5);
                voxels_ii->SetLineColor(linecolor);
                //voxels_ii->SetLineColorAlpha(kMagenta, 0.5);
            }

            
            if (padtype != _Polar)
            {
                x_node = (r_start + r_end) / 2.;
                y_node = (phi_start + phi_end) / 2. / TMath::RadToDeg();
            }
            else
            {
                x_node = 0.;
                y_node = 0.;
            }
            z_node = m_enable3D ? zidx0 : 0;
            //auto node_ii = new TNode(Form("node_%d", ii), "node", Form("voxels_%d", ii), x_node, y_node, zidx0);
            auto node_ii = new TNode(Form("node_%d", ii), "node", Form("voxels_%d", ii), x_node, y_node, z_node);
            vzbinSig.push_back(zidx0);
            //std::printf("[debug] PixelIdx=%d, ZbinIdx=%d, x_node=%.4f, y_node=%.4f\n", xyidx0, zidx0, x_node, y_node);
        }
    }

    int Counter = 0;
    auto min_ele = std::min_element(vzbinSig.begin(), vzbinSig.end());
    auto max_ele = std::max_element(vzbinSig.begin(), vzbinSig.end());

    // bkg
    for (long long i_entry = 0; i_entry < m_tbkg->GetEntries(); ++i_entry)
    {
        m_tbkg->GetEntry(i_entry);
        if (m_zbinIdx >= (*min_ele) && m_zbinIdx <= (*max_ele))
        {
            int bkgxyidx = static_cast<int>(m_pixelIdx);
            auto padvertex_ii = helper->GetPadVertex(bkgxyidx);
            double r_start = std::get<0>(padvertex_ii);
            double r_end = std::get<1>(padvertex_ii);
            double phi_start = std::get<2>(padvertex_ii) * TMath::RadToDeg();
            double phi_end = std::get<3>(padvertex_ii) * TMath::RadToDeg();
            bool conditionPhi = false;
            if (padtype == _Polar)
                conditionPhi = (phi_start >= (trkPhi - m_phiRange)) && (phi_end <= (trkPhi + m_phiRange));
            else
            {
                double phip = TMath::ATan2(phi_start / TMath::RadToDeg(), r_start);

                while (phip < 0) phip += TMath::TwoPi();
                phip *= TMath::RadToDeg();
                conditionPhi = (phip >= (trkPhi - m_phiRange)) && (phip <= (trkPhi + m_phiRange));
            }

            if (conditionPhi)
            {
                if (padtype == _Polar)
                {
                    auto voxels_ii = new TTUBS(Form("bkgvoxels_%lld", i_entry), "void", "vacuum", r_start, r_end, 0.5, phi_start, phi_end);
                    voxels_ii->SetLineColor(kGreen);
                }
                else
                {
                    auto voxels_ii = new TBRIK(Form("bkgvoxels_%lld", i_entry), "void", "vacuum", helper->GetPadWidth() * 0.5, helper->GetPadWidth() * 0.5, 0.5);
                    voxels_ii->SetLineColor(kGreen);
                }

                if (padtype != _Polar)
                {
                    x_node = (r_start + r_end) / 2.;
                    y_node = (phi_start + phi_end) / 2. / TMath::RadToDeg();
                }
                else
                {
                    x_node = 0.;
                    y_node = 0.;
                }
                z_node = m_enable3D ? m_zbinIdx : 0.;
                auto node_ii = new TNode(Form("bkgnode_%lld", i_entry), "node", Form("bkgvoxels_%lld", i_entry), x_node, y_node, z_node);
                Counter++;
            }
        }
    }

    std::printf("[debug]===> Counter=%d, min=%d, max=%d\n", Counter, *min_ele, *max_ele);
    topNode->SetVisibility(0);

    return topNode;
}

void PlotVoxelsGL::AddAnnotation(const char* txt)
{
    TGLAnnotation *ann = new TGLAnnotation(m_glviewer, txt, 0.05, 0.92);
    ann->SetTextSize(0.05);
    m_glviewer->RefreshPadEditor(m_glviewer);
    m_glviewer->DoDraw();
}
