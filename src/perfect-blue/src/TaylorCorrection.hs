{-# LANGUAGE FlexibleContexts #-}

-- Taylor Series Correction System with Information Criterion
-- Based on finite-difference Taylor expansions and least-squares correction

module TaylorCorrection where

import Data.List (transpose)

-- Types
type Buffer = [Double]           -- Buffered sequence points
type TaylorCoeffs = [Double]     -- Taylor series coefficients
type Matrix = [[Double]]         -- 2D matrix
type Vector = [Double]           -- 1D vector

-- Configuration
data ModelConfig = ModelConfig
    { historyLength :: Int       -- H: full history buffer size
    , lowerOrder :: Int          -- K: reduced Taylor order (K < H)
    , cvLength :: Int            -- Cross-validation segment size
    , deltaT :: Double           -- Time step spacing
    , lambda :: Double           -- Complexity penalty weight
    } deriving Show

-- ============================================================================
-- Matrix Operations
-- ============================================================================

-- Matrix multiplication
matMul :: Matrix -> Matrix -> Matrix
matMul a b = [[sum $ zipWith (*) row col | col <- transpose b] | row <- a]

-- Matrix-vector multiplication
matVec :: Matrix -> Vector -> Vector
matVec m v = [sum $ zipWith (*) row v | row <- m]

-- Vector-matrix multiplication
vecMat :: Vector -> Matrix -> Vector
vecMat v m = [sum $ zipWith (*) v col | col <- transpose m]

-- Matrix transpose
matTranspose :: Matrix -> Matrix
matTranspose = transpose

-- Matrix inversion using Gauss-Jordan (from previous code)
matInvert :: Matrix -> Maybe Matrix
matInvert m
    | rows /= cols = Nothing
    | otherwise = gaussJordan augmented
  where
    rows = length m
    cols = length (head m)
    identity = [[ if i == j then 1 else 0 | j <- [0..cols-1]] | i <- [0..rows-1]]
    augmented = zipWith (++) m identity

gaussJordan :: Matrix -> Maybe Matrix
gaussJordan m = do
    reduced <- forwardElimination m 0
    result <- backSubstitution reduced (length reduced - 1)
    return $ map (drop (length result)) result

forwardElimination :: Matrix -> Int -> Maybe Matrix
forwardElimination m row
    | row >= length m = Just m
    | otherwise = do
        m' <- pivotRow m row
        let m'' = eliminateColumn m' row
        forwardElimination m'' (row + 1)

pivotRow :: Matrix -> Int -> Maybe Matrix
pivotRow m row =
    case findPivot m row row of
        Nothing -> Nothing
        Just pivotIdx -> Just $ swapRows m row pivotIdx

findPivot :: Matrix -> Int -> Int -> Maybe Int
findPivot m col row
    | row >= length m = Nothing
    | abs (m !! row !! col) > 1e-10 = Just row
    | otherwise = findPivot m col (row + 1)

swapRows :: Matrix -> Int -> Int -> Matrix
swapRows m i j =
    [ if k == i then m !! j
      else if k == j then m !! i
      else m !! k
    | k <- [0..length m - 1]]

eliminateColumn :: Matrix -> Int -> Matrix
eliminateColumn m row =
    [ if i <= row then m !! i
      else let pivot = m !! row !! row
               factor = (m !! i !! row) / pivot
           in zipWith (\a b -> a - factor * b) (m !! i) (m !! row)
    | i <- [0..length m - 1]]

backSubstitution :: Matrix -> Int -> Maybe Matrix
backSubstitution m row
    | row < 0 = Just m
    | abs (m !! row !! row) < 1e-10 = Nothing
    | otherwise = do
        let m' = scaleRow m row
            m'' = eliminateAbove m' row
        backSubstitution m'' (row - 1)

scaleRow :: Matrix -> Int -> Matrix
scaleRow m row =
    [ if i /= row then m !! i
      else let pivot = m !! row !! row
           in map (/ pivot) (m !! i)
    | i <- [0..length m - 1]]

eliminateAbove :: Matrix -> Int -> Matrix
eliminateAbove m row =
    [ if i >= row then m !! i
      else let factor = m !! i !! row
           in zipWith (\a b -> a - factor * b) (m !! i) (m !! row)
    | i <- [0..length m - 1]]

-- ============================================================================
-- Vandermonde Matrix Construction
-- ============================================================================

factorial :: Int -> Double
factorial 0 = 1
factorial n = fromIntegral n * factorial (n - 1)

-- Construct Vandermonde design matrix for Taylor series
-- Rows correspond to buffer indices, columns to Taylor terms
vandermonde :: [Int] -> Double -> Int -> Matrix
vandermonde indices dt order =
    [[ (fromIntegral idx * dt)^i / factorial i 
     | i <- [0..order]] 
    | idx <- indices]

-- ============================================================================
-- Taylor Coefficient Computation
-- ============================================================================

-- Compute Taylor coefficients from buffer using least squares
-- This solves: minimize ||M*t - buffer||^2
computeTaylorCoeffs :: Buffer -> Double -> Int -> Maybe TaylorCoeffs
computeTaylorCoeffs buffer dt order = do
    let n = length buffer
        indices = [0..n-1]
        m = vandermonde indices dt order
        mt = matTranspose m
        mtm = matMul mt m
    
    mtmInv <- matInvert mtm
    let mtb = matVec mt buffer
        coeffs = matVec mtmInv mtb
    return coeffs

-- ============================================================================
-- Correction Matrix Computation
-- ============================================================================

-- Compute the correction matrix B that maps full-order to lower-order coefficients
-- B = (M_K^T M_K)^(-1) M_K^T M_H
-- This is precomputed and data-independent
computeCorrectionMatrix :: ModelConfig -> IO (Maybe Matrix)
computeCorrectionMatrix config = do
    -- CV indices are the time steps for the cross-validation segment
    -- Starting from step after history ends
    let histLen = historyLength config + 1  -- number of history points
        cvLen = cvLength config
        cvIndices = [histLen .. histLen + cvLen - 1]
        dt = deltaT config
        fullOrder = historyLength config
        lowOrder = lowerOrder config
    
    putStrLn $ "CV indices: " ++ show cvIndices
    putStrLn $ "Full order: " ++ show fullOrder ++ ", Low order: " ++ show lowOrder
    
    -- Design matrices for CV segment
    let mH = vandermonde cvIndices dt fullOrder
        mK = vandermonde cvIndices dt lowOrder
    
    putStrLn $ "M_H dimensions: " ++ show (length mH) ++ "x" ++ show (length $ head mH)
    putStrLn $ "M_K dimensions: " ++ show (length mK) ++ "x" ++ show (length $ head mK)
    
    putStrLn "M_K:"
    mapM_ print mK
    
    -- Compute (M_K^T M_K)^(-1)
    let mKt = matTranspose mK
        mKtmK = matMul mKt mK
    
    putStrLn "\nM_K^T M_K:"
    mapM_ print mKtmK
    
    case matInvert mKtmK of
        Nothing -> do
            putStrLn "Failed to invert M_K^T M_K"
            return Nothing
        Just mKtmKinv -> do
            putStrLn "\n(M_K^T M_K)^(-1):"
            mapM_ print mKtmKinv
            
            -- Compute B = (M_K^T M_K)^(-1) M_K^T M_H
            let mKtmH = matMul mKt mH
                b = matMul mKtmKinv mKtmH
            return $ Just b

-- ============================================================================
-- Prediction and Correction
-- ============================================================================

-- Apply correction to lower-order coefficients
-- T_K^corrected = T_K^hist + B * (T_H^full - [T_K^hist; 0...0])
applyCorrection :: TaylorCoeffs -> TaylorCoeffs -> Matrix -> TaylorCoeffs
applyCorrection tLowHist tFull bMatrix = 
    zipWith (+) tLowHist correction
  where
    -- Pad lower-order coeffs with zeros to match full length
    tLowPadded = tLowHist ++ replicate (length tFull - length tLowHist) 0
    -- Compute residual
    residual = zipWith (-) tFull tLowPadded
    -- Apply correction matrix
    correction = matVec bMatrix residual

-- Predict future values using Taylor coefficients
-- f(t) = sum_{i=0}^k c_i * (dt*n)^i / i!
predict :: TaylorCoeffs -> Double -> Int -> [Double]
predict coeffs dt nSteps =
    [ sum [ c * ((fromIntegral step * dt)^i / factorial i) 
          | (i, c) <- zip [0..] coeffs ]
    | step <- [1..nSteps] ]

-- ============================================================================
-- Cross-Validation and Information Criterion
-- ============================================================================

-- Compute prediction error on validation data
predictionError :: Buffer -> TaylorCoeffs -> Double -> Double
predictionError validation coeffs dt =
    sum $ zipWith (\actual pred -> (actual - pred)^2) validation predictions
  where
    nSteps = length validation
    predictions = predict coeffs dt nSteps

-- Compute information criterion: IC = E_cv + λ * complexity
informationCriterion :: Double -> Int -> Double -> Double
informationCriterion errorCV order lambda =
    errorCV + lambda * fromIntegral order

-- ============================================================================
-- Complete Workflow
-- ============================================================================

-- Process a single buffer sequence with correction
processBuffer :: ModelConfig -> Buffer -> IO (Maybe (TaylorCoeffs, TaylorCoeffs, Double))
processBuffer config buffer = do
    let histLen = historyLength config
        lowOrder = lowerOrder config
        dt = deltaT config
        
        -- Split into history and CV
        history = take (histLen + 1) buffer
        cv = drop (histLen + 1) buffer
    
    putStrLn $ "History: " ++ show history
    putStrLn $ "CV: " ++ show cv
    
    -- Compute full-order Taylor on history
    case computeTaylorCoeffs history dt histLen of
        Nothing -> do
            putStrLn "Failed to compute full-order Taylor coefficients"
            return Nothing
        Just tFull -> do
            putStrLn $ "Full Taylor coeffs: " ++ show tFull
            
            -- Compute lower-order Taylor on history
            case computeTaylorCoeffs history dt lowOrder of
                Nothing -> do
                    putStrLn "Failed to compute lower-order Taylor coefficients"
                    return Nothing
                Just tLowHist -> do
                    putStrLn $ "Lower Taylor coeffs: " ++ show tLowHist
                    
                    -- Get correction matrix
                    bMatrixResult <- computeCorrectionMatrix config
                    case bMatrixResult of
                        Nothing -> do
                            putStrLn "Failed to compute correction matrix"
                            return Nothing
                        Just bMatrix -> do
                            putStrLn "\nCorrection matrix B:"
                            mapM_ print bMatrix
                            
                            -- Apply correction
                            let tLowCorrected = applyCorrection tLowHist tFull bMatrix
                            
                            -- Compute CV error
                            let errorCV = predictionError cv tLowCorrected dt
                                ic = informationCriterion errorCV lowOrder (lambda config)
                            
                            return $ Just (tLowHist, tLowCorrected, ic)

-- ============================================================================
-- Example Usage
-- ============================================================================

exampleConfig :: ModelConfig
exampleConfig = ModelConfig
    { historyLength = 3
    , lowerOrder = 2
    , cvLength = 2
    , deltaT = 1.0
    , lambda = 1.0
    }

-- Generate test data: simple polynomial
testData :: [Double]
testData = [x^2 + 2*x + 1 | x <- [0..6]]

main :: IO ()
main = do
    putStrLn "Taylor Series Correction System"
    putStrLn "================================\n"
    
    putStrLn "Configuration:"
    print exampleConfig
    putStrLn ""
    
    putStrLn "Test data:"
    print testData
    putStrLn ""
    
    result <- processBuffer exampleConfig testData
    case result of
        Nothing -> putStrLn "Error: Unable to process buffer"
        Just (tLowHist, tLowCorrected, ic) -> do
            putStrLn "\nLower-order coefficients (history only):"
            print tLowHist
            putStrLn ""
            
            putStrLn "Corrected lower-order coefficients:"
            print tLowCorrected
            putStrLn ""
            
            putStrLn $ "Information Criterion: " ++ show ic
            putStrLn ""
            
            -- Make predictions
            let predictions = predict tLowCorrected (deltaT exampleConfig) 2
            putStrLn "Predictions for next 2 steps:"
            print predictions
            putStrLn ""
            
            let actual = drop 5 testData
            putStrLn "Actual values:"
            print actual