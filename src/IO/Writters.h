#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include <Eigen/Dense>

namespace IO
{

// template <bool un_signed = true
template <typename PointCloudDiff>
void savePointCloud( PointCloudDiff& pointCloud, std::string & filename )
{
    // if (un_signed){
    //   switch_to_abs ();      
    // }

    std::ofstream ofs( filename, std::ofstream::out );
    ofs << "# x y z Gauss_Curvature Mean_Curvature nx ny nz k1 k2 d1x d1y d1z d2x d2y d2z\n";
    std::cout << "Number of points : " << pointCloud.points.size() << std::endl;
    for ( auto i = 0u; i < pointCloud.points.size(); ++i )
    {
        ofs << pointCloud.points[ i ][ 0 ] << " " << pointCloud.points[ i ][ 1 ] << " " << pointCloud.points[ i ][ 2 ];
        ofs << " " << pointCloud.gauss[ i ] << " " << pointCloud.mean[ i ] << " ";
        ofs << pointCloud.normals[ i ][ 0 ] << " " << pointCloud.normals[ i ][ 1 ] << " " << pointCloud.normals[ i ][ 2 ] << " ";
        ofs << pointCloud.k1[ i ] << " " << pointCloud.k2[ i ] << " ";
        ofs << pointCloud.v1[ i ][ 0 ] << " " << pointCloud.v1[ i ][ 1 ] << " " << pointCloud.v1[ i ][ 2 ] << " ";
        ofs << pointCloud.v2[ i ][ 0 ] << " " << pointCloud.v2[ i ][ 1 ] << " " << pointCloud.v2[ i ][ 2 ];
        ofs << std::endl;
    }
    ofs.close();
}

// template <bool un_signed = true
template <typename PointCloudDiff>
void savePointCloudAsPositions( PointCloudDiff& pointCloud, std::string & filename )
{
    // if (un_signed){
    //   switch_to_abs ();      
    // }

    std::ofstream ofs( filename, std::ofstream::out );
    ofs << "# x y z\n";
    std::cout << "Number of points : " << pointCloud.points.size() << std::endl;
    for ( auto i = 0u; i < pointCloud.points.size(); ++i )
    {
        ofs << pointCloud.points[ i ][ 0 ] << " " << pointCloud.points[ i ][ 1 ] << " " << pointCloud.points[ i ][ 2 ];
        ofs << std::endl;
    }
    ofs.close();
}

// template <bool un_signed = true>
template <typename PointCloudDiff, typename _Scalar>
void savePointCloudAsErrors( PointCloudDiff& pointCloud, std::string & filename )
{
    // if (un_signed){
    //   switch_to_abs ();
    //   other.switch_to_abs ();      
    // }

    std::ofstream ofs( filename, std::ofstream::out );
    ofs << "# x y z iShape Gauss_Curvature Mean_Curvature k1 k2 d1 d2 normal pos\n";
    for ( auto i = 0u; i < pointCloud.points.size(); ++i )
    {

        _Scalar errPos  = pointCloud.errorPos[ i ];
        _Scalar errNorm = pointCloud.errorNormal[ i ];
        _Scalar errMean = pointCloud.errorMean[ i ];
        _Scalar errGauss = pointCloud.errorGauss[ i ];
        _Scalar errK1 = pointCloud.errorK1[ i ];
        _Scalar errK2 = pointCloud.errorK2[ i ];
        _Scalar errD1 = pointCloud.errorD1[ i ];
        _Scalar errD2 = pointCloud.errorD2[ i ];
        _Scalar errSIndex = pointCloud.errorShapeIndex[ i ];

        ofs << pointCloud.points[ i ][ 0 ] << " " << pointCloud.points[ i ][ 1 ] << " " << pointCloud.points[ i ][ 2 ];
        ofs << " " << errSIndex << " " << errGauss << " " << errMean << " ";
        ofs << errK1 << " " << errK2 << " ";
        ofs << errD1 << " " << errD2 << " ";
        ofs << errNorm << " ";
        ofs << errPos;
        ofs << "\n";
    }
    ofs.close();
}

template <typename PointCloudDiff, typename _Scalar>
std::string getStats (const PointCloudDiff & estimation, std::string & name, _Scalar radius, _Scalar noise_position, _Scalar noise_normal, _Scalar flip_ratio, bool absolute_err = true ) {

    _Scalar stable_ratio = 0;
    
    for (auto & stable : estimation.non_stable_idx ) { 
        stable_ratio += stable;
    }
    stable_ratio /= estimation.points.size();

    std::cout << "# nbPoints radius noise-position noise-normal flip-normal {avg,max,variance} for [nbNeighbors,mean,gauss,k1,k2,d1,d2,pos,idxShape,normal,timings] non_stable_ratio and variant-name\n";
    std::string variant = name;

    std::vector <_Scalar> errors = { estimation.errorMean.mean(), estimation.errorMean.max(), estimation.errorMean.unbiasedVariance(),
                                     estimation.errorGauss.mean(), estimation.errorGauss.max(), estimation.errorGauss.unbiasedVariance(),
                                     estimation.errorK1.mean(), estimation.errorK1.max(), estimation.errorK1.unbiasedVariance(),
                                     estimation.errorK2.mean(), estimation.errorK2.max(), estimation.errorK2.unbiasedVariance(),
                                     estimation.errorD1.mean(), estimation.errorD1.max(), estimation.errorD1.unbiasedVariance(),
                                     estimation.errorD2.mean(), estimation.errorD2.max(), estimation.errorD2.unbiasedVariance(),
                                     estimation.errorPos.mean(), estimation.errorPos.max(), estimation.errorPos.unbiasedVariance(),
                                     estimation.errorShapeIndex.mean(), estimation.errorShapeIndex.max(), estimation.errorShapeIndex.unbiasedVariance(), 
                                     estimation.errorNormal.mean(), estimation.errorNormal.max(), estimation.errorNormal.unbiasedVariance() } ;

    if ( ! absolute_err ) { // If not absolute error, then the errors are MSE, so compute the square root
        for (auto & e : errors) {
            e = std::sqrt(e);
        }
    }
        

    std::ostringstream stats;
    stats        << estimation.points.size() << " " << radius << " " << noise_position << " " << noise_normal << " " << flip_ratio << " "
                 << estimation.statNeighbors.mean() << " " << estimation.statNeighbors.max() << " " << estimation.statNeighbors.unbiasedVariance() << " " ;
    
    for (auto & e : errors) {
        stats    << e << " ";
    }
    stats        << estimation.statTimings.mean() << " " << estimation.statTimings.max() << " " << estimation.statTimings.unbiasedVariance() << " "
                 << stable_ratio << " " << variant << std::endl;
    return stats.str();
}

// template <typename PointCloudDiff>
// void saveStats( PointCloudDiff& pointCloud, std::string & stats, std::string & filename )
// {
//     std::ofstream ofs( statsFilename, std::ios::out | std::ios::app );
//     ofs << stats;
//     ofs.close();
// }

} // namespace IO
