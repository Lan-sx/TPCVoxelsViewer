/*********************************************************************
 * Author           : Lan-sx
 * Email            : shexin@ihep.ac.cn
 * Created          : 2025-03-31 16:58
 * Filename         : TPCVoxelHelper.h
 * Description      : A helper class for TPC Voxels analysis
 * ******************************************************************/
#ifndef __TPCVOXELHELPER__
#define __TPCVOXELHELPER__ 1
//std
#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <tuple>
#include <numeric>

//ROOT CERN
#include "TMath.h"

//const expr
constexpr double _TPCinnerR = 63.5;   // [cm]
constexpr double _TPCouterR = 175.0;  // [cm]
constexpr double _TPCMaxDftL = 276.0; // [cm]
constexpr double _Boxsize = 180.*2;   // [cm]

//Cartesian, default, square pixel layout
struct TPCVoxels
{
    int PixelIdx;
    int ZbinIdx;

    TPCVoxels(int pixelIndex=0, int zbinIndex=0) : PixelIdx(pixelIndex), ZbinIdx(zbinIndex) {}

    bool operator<(const TPCVoxels& others) const
    {
        if(others.PixelIdx == PixelIdx)
        {
            return ZbinIdx < others.ZbinIdx;
        }
        else
            return PixelIdx < others.PixelIdx;
    }

};

//Polar Row by Row pixel/pad layout 
//struct TPCVoxelsRPhi
//{
//    int RowIdx;
//    int AngleIdx;
//    int ZbinIdx;
//
//    TPCVoxelsRPhi(int rowIndex=0, int angleIndex=0, int zbinIndex=0) : RowIdx(rowIndex), AngleIdx(angleIndex), ZbinIdx(zbinIndex) {}
//
//    bool operator<(const TPCVoxelsRPhi& others) const 
//    {
//        if(ZbinIdx < others.ZbinIdx) return true;
//        if(ZbinIdx > others.ZbinIdx) return false;
//
//        if(RowIdx < others.RowIdx) return true;
//        if(RowIdx > others.RowIdx) return false;
//
//        return AngleIdx < others.AngleIdx;
//    }
//};

//enum of Voxel type
enum VoxelType { _Cartesian=0, _Polar};

class TPCVoxelHelper
{
public:
    //methods
    static TPCVoxelHelper* GetInstance()
    {
        if(m_Instance == nullptr)
        {
            m_Instance = new TPCVoxelHelper();
        }
        return m_Instance;
    }

    void InitialVoxel(int voxeltype, double padlength, double padwidth);

    void Destruct()
    {
        if(m_Instance != nullptr)
        {
            delete m_Instance;
            m_Instance = nullptr;
        }
    }
    
    //Setter Methods
    inline void SetTPCVoxel(int type, double length, double width);  

    //Getter Methods
    int GetVoxelType() { return m_voxelType; }
    double GetPadLenth() { return m_padLength; }
    double GetPadWidth() { return m_padWidth; } 
    int GetNrowsInRadius() { return m_NpadperRow.size(); }
    const std::vector<int> &GetNpadperRowVector() const{ return m_NpadperRow; }
    const int GetTotalpadsInxyPlane(); 
    int GetGlobalxyIdx(double xx, double yy);
    //===========================================================
    // for simple square pad layout, padindex1 -> row idx, padindex2 -> col idx
    // for polar pad layout, padindex1 -> row idx, padindex2 -> phi idx
    int GetGlobalxyIdxFromPadIndex12(int padindex1, int padindex2);
    
    //===========================================================
    //Get Pad index at x-y plane from globalxyIdx
    //for simple square pad layout, return row col index
    //for polar pad layout, return row and phi index
    std::pair<int, int> GetPadIndex12(int globalxyIdx);
    //===========================================================
    //Get Pad vertex at x-y plane from globalxyIdx
    //for simple square pad layout, return pad center
    //for polar pad layout, return r1 r2, phi1,phi2
    std::tuple<double, double, double, double> GetPadVertex(int globalxyIdx);

protected:
    int GetGlobalxyIdx_cartesian(double xx, double yy);
    int GetGlobalxyIdx_polar(double xx, double yy);
    
    // used for Row by Row pad layout
    void CreateNpadperRowVec();
    

private:
    TPCVoxelHelper() : m_voxelType(_Cartesian), m_padLength(0.05), m_padWidth(0.05) {};
    static TPCVoxelHelper* m_Instance;

private:
    int m_voxelType;
    double m_padLength, m_padWidth;
    std::vector<int> m_NpadperRow;
};

inline void TPCVoxelHelper::SetTPCVoxel(int type, double length, double width)
{
    this->InitialVoxel(type,length,width);
}

#endif
