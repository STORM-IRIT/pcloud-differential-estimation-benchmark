# Ponca estimators

You can simply add ponca concept into the directory `classes` (if it doesn't exist yet) and use it in the `estimators.h` file.

You just need to add the function into the `estimators.h` like this example :

```cpp
DifferentialQuantities<Scalar> estimateDifferentialQuantitiesWith_YOUR_METHOD_() {
    using basket_YOUR_METHOD_Fit = Ponca::Basket<PPAdapter, SmoothWeightFunc, Ponca::_YOUR_METHOD_Fit>;
    using basketDiff__YOUR_METHOD_Fit = Ponca::BasketDiff<
            basket__YOUR_METHOD_Fit,
            Ponca::DiffType::FitSpaceDer,
            Ponca::_YOUR_METHOD_Der,
            Ponca::CurvatureEstimatorBase, Ponca::NormalDerivativesCurvatureEstimator>;
    return estimateDifferentialQuantities_impl<basketDiff__YOUR_METHOD_Fit>("_YOUR_METHOD_");
}
```

With that, you need to manualy add this new feature into the `compute_ponca_estimators.cpp` file.

## 3D Ellipsoid estimator

This estimator computes the algebraic function by directly fitting the ellipsoid using the information about the normals.

## Explanation about Cylinder / 2D Ellipsoid estimators

These estimators are base on the same pipeline comming from the [PC-MLS][1] (Parabolic-Cylindrical Moving Least Squares Surfaces).
 
### Cylindrical shape

- The Base Cylinder estimator compute a covariance plane and fit its parabolic-cylindrical like the paper. 
- The Fully Oriented Cylinder estimator compute a mean plane and directly fit its parabolic-cylindrical by using the information about the normals.
- The Near Oriented Cylinder estimator compute a Covariance plane and directly fit its parabolic-cylindrical by using the information about the normals.
- The Base Oriented Cylinder estimator compute a mean plane and fit its parabolic-cylindrical like the paper.

### 2D Ellipsoid shape

If we take into account the paper, we fit a plane, then, we fit the height value to get a 2D Ellipsoid, and after, we force this shape to not quadraticly vary in a certain direction to obtain a cylindrical shape. 

Here, we stop the computation after 2D ellipsoid fitting.

- The Base Oriented 2D Ellipsoid estimator compute a mean plane and fit its 2D Ellipsoid like the paper.
- The Fully Oriented 2D Ellipsoid estimator compute a mean plane and directly fit its 2D Ellipsoid by using the information about the normals.



[1]:https://inria.hal.science/hal-01169572/file/main.pdf