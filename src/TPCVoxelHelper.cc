/*********************************************************************
 * Author           : Lan-sx
 * Email            : shexin@ihep.ac.cn
 * Created          : 2025-04-01 10:29
 * Filename         : TPCVoxelHelper.cc
 * Description      : 
 * ******************************************************************/
#include "TPCVoxelHelper.hh"

TPCVoxelHelper* TPCVoxelHelper::m_Instance = nullptr;  // definition

void TPCVoxelHelper::InitialVoxel(int voxeltype, double padlength, double padwidth)
{
    if (m_Instance)
    {
        m_padLength = padlength;
        m_padWidth = padwidth;
        switch (voxeltype)
        {
        case _Cartesian:
            m_voxelType = _Cartesian;
            break;
        case _Polar:
            m_voxelType = _Polar;
            m_NpadperRow.clear();
            this->CreateNpadperRowVec();
            break;
        default:
            std::printf("[warning] unknown type!!!\n");
            break;
        }
    }
}

int TPCVoxelHelper::GetGlobalxyIdx(double xx, double yy)
{
    switch (m_voxelType)
    {
    case _Cartesian:
        return this->GetGlobalxyIdx_cartesian(xx,yy);
        break;
    case _Polar:
        return this->GetGlobalxyIdx_polar(xx,yy);
        break;
    default:
        return -1;
        break;
    }
}

int TPCVoxelHelper::GetGlobalxyIdxFromPadIndex12(int padindex1, int padindex2)
{
    if(padindex1 < 0 || padindex2 < 0)
        return -1;
    else
    {
        if (m_voxelType == _Polar)
        {
            //check row idx range
            if(padindex1 >= this->GetNrowsInRadius())
                return -1;
            //check phi idx range, phi idx must be small than the Npad in padindex1-th row
            if(padindex2 >= m_NpadperRow.at(padindex1))
                return -1;
            
            auto sumRow = std::accumulate(m_NpadperRow.begin(), m_NpadperRow.begin() + padindex1, 0.);
            int globalIdx = padindex2 + static_cast<int>(sumRow);
            return globalIdx;
        }
        else if (m_voxelType == _Cartesian)
        {
            int nBox = static_cast<int>(_Boxsize / m_padWidth);
            if (padindex1 >= nBox || padindex2 >= nBox)
                return -1;
            else
                return (padindex1 * nBox + padindex2);
        }
        else
            return -1;
    }
}

int TPCVoxelHelper::GetGlobalxyIdx_cartesian(double xx, double yy)
{
    double radius = TMath::Sqrt(xx*xx+yy*yy);
    if( radius < _TPCinnerR || radius > _TPCouterR)
        return -1;
    int nBox = static_cast<int>(_Boxsize/m_padWidth);
    int row_id = static_cast<int>(std::floor(yy/m_padWidth)+nBox/2);
    int col_id = static_cast<int>(std::floor(xx/m_padWidth)+nBox/2);
    //std::printf("[debug] nBox=%d, %d %d %.2f padwidth=%.2f\n", nBox, row_id, col_id,std::floor(xx/m_padWidth),m_padWidth);
    row_id =  (row_id >= nBox) ? -1 : row_id;
    col_id =  (col_id >= nBox) ? -1 : col_id;

    if(row_id >= 0 && col_id >= 0)
        return row_id*nBox + col_id;
    else 
        return -1;
}

int TPCVoxelHelper::GetGlobalxyIdx_polar(double xx, double yy)
{
    double radius = TMath::Sqrt(xx*xx+yy*yy);
    if( radius < _TPCinnerR || radius > _TPCouterR)
        return -1;
    
    double phi = TMath::ATan2(yy,xx);
    while(phi < 0.) phi += 2*TMath::Pi();
    
    //check row_id range
    int row_id = static_cast<int>((radius-_TPCinnerR)/m_padLength); 
    if(row_id >= this->GetNrowsInRadius())
        return -1;
    
    auto radius_i = _TPCinnerR + row_id*m_padLength;

    double dphi = 2*TMath::Pi()/m_NpadperRow.at(row_id);
    int phi_id = static_cast<int>(phi/dphi);

    auto sumRow = std::accumulate(m_NpadperRow.begin(), m_NpadperRow.begin()+row_id,0.);
    int globalIdx = phi_id + static_cast<int>(sumRow); 

    return globalIdx;
}

void TPCVoxelHelper::CreateNpadperRowVec()
{
    int Nrow = static_cast<int>((_TPCouterR-_TPCinnerR)/m_padLength);
    m_NpadperRow.resize(Nrow);

    for(int row_i=0; row_i < Nrow; ++row_i)
    {
        double radius_i = _TPCinnerR + row_i*m_padLength;
        auto Npad = static_cast<int>(2*TMath::Pi()*radius_i/m_padWidth);
        m_NpadperRow.at(row_i) = Npad;
    }
}

std::pair<int, int> TPCVoxelHelper::GetPadIndex12(int globalxyIdx)
{
    if (m_voxelType == _Polar)
    {
        int globalCounter1 = 0, globalCounter2 = 0;

        int rowindex = -1;
        int phiindex = -1;

        for (size_t row_i = 0; row_i < m_NpadperRow.size(); ++row_i)
        {
            globalCounter1 += m_NpadperRow.at(row_i);
            if (globalxyIdx >= globalCounter2 && globalxyIdx <= globalCounter1)
            {
                rowindex = row_i;
                phiindex = globalxyIdx - globalCounter2;
            }
            globalCounter2 = globalCounter1;
        }

        return std::make_pair(rowindex, phiindex);
    }
    else if (m_voxelType == _Cartesian)
    {
        int nBox = static_cast<int>(_Boxsize/m_padWidth);
        int rowindex = globalxyIdx/nBox;
        int colindex = globalxyIdx%nBox;

        return std::make_pair(rowindex,colindex); 
    }
    else
    {
        std::printf("[warning] unknown pad layout type!\n");
        return std::make_pair(-1,-1);
    }
}

std::tuple<double, double, double, double> TPCVoxelHelper::GetPadVertex(int globalxyIdx)
{
    auto padidx12 = this->GetPadIndex12(globalxyIdx);

    if (m_voxelType == _Polar)
    {
        if (padidx12.first >= 0 && padidx12.second >= 0)
        {
            double r_start = _TPCinnerR + padidx12.first * m_padLength;
            double r_end = r_start + m_padLength;

            double dphi = TMath::TwoPi() / m_NpadperRow.at(padidx12.first);
            double angle_start = padidx12.second * dphi;
            double angle_end = angle_start + dphi;
            return std::make_tuple(r_start, r_end, angle_start, angle_end);
        }
        else
            return std::make_tuple(0.,0.,0.,0.);
    }
    else if (m_voxelType == _Cartesian)
    {
        int nBox = static_cast<int>(_Boxsize/m_padWidth);
        double x_low = (-nBox/2 + padidx12.second)*m_padWidth;
        double x_up = x_low + m_padWidth;
        double y_low = (-nBox/2 + padidx12.first)*m_padWidth;
        double y_up = y_low + m_padWidth;

        return std::make_tuple(x_low,x_up,y_low,y_up); 
    }else 
    {
        return std::make_tuple(0., 0., 0., 0.);
    }
}

const int TPCVoxelHelper::GetTotalpadsInxyPlane()
{
    if(m_voxelType == _Cartesian)
    {
        int TotalPads = static_cast<int>(TMath::Pi()*(_TPCouterR*_TPCouterR-_TPCinnerR*_TPCinnerR)/(m_padWidth*m_padWidth));
        return TotalPads;
    }
    else if(m_voxelType == _Polar)
    {
        int TotalPads = std::accumulate(m_NpadperRow.begin(), m_NpadperRow.end(), 0.);
        return TotalPads;
    }
    else
    {
        return 1;
    }

}
