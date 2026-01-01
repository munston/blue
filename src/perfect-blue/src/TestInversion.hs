module TestInversion where

import MatrixInversion

main :: IO ()
main = do
    putStrLn "=== Matrix Inversion Tests ==="
    putStrLn ""
    
    -- Test 1: Invertible 3x3 matrix
    putStrLn "Test 1: Invertible 3x3 matrix"
    putStrLn "------------------------------"
    let m1 = [[4, 7, 2],
              [1, 3, 1],
              [2, 5, 3]] :: Matrix
    
    putStrLn "Original matrix:"
    printMatrix m1
    
    case invertMatrix m1 of
        Just minv -> do
            putStrLn "\nInverse matrix:"
            printMatrix minv
            
            -- Verify: M * M^-1 should be identity
            let identity = cleanZeros $ matMul m1 minv
            putStrLn "\nVerification (M * M^-1):"
            printMatrix identity
            
        Nothing -> putStrLn "Matrix is not invertible!"
    
    putStrLn ""
    
    -- Test 2: Singular matrix
    putStrLn "Test 2: Singular matrix"
    putStrLn "------------------------"
    let m2 = [[1, 2, 3],
              [2, 4, 6],
              [1, 2, 3]] :: Matrix
    
    putStrLn "Singular matrix:"
    printMatrix m2
    
    case invertMatrix m2 of
        Just _ -> putStrLn "\nUnexpectedly got an inverse"
        Nothing -> putStrLn "\nMatrix is not invertible (as expected)"
    
    putStrLn ""
    
    -- Test 3: 2x2 matrix
    putStrLn "Test 3: Simple 2x2 matrix"
    putStrLn "--------------------------"
    let m3 = [[4, 7],
              [2, 6]] :: Matrix
    
    putStrLn "Original matrix:"
    printMatrix m3
    
    case invertMatrix m3 of
        Just minv -> do
            putStrLn "\nInverse matrix:"
            printMatrix minv
            
            let identity = cleanZeros $ matMul m3 minv
            putStrLn "\nVerification (M * M^-1):"
            printMatrix identity
            
        Nothing -> putStrLn "Matrix is not invertible!"
    
    putStrLn ""
    
    -- Test 4: Identity matrix
    putStrLn "Test 4: Identity matrix"
    putStrLn "-----------------------"
    let m4 = [[1, 0, 0],
              [0, 1, 0],
              [0, 0, 1]] :: Matrix
    
    putStrLn "Original matrix:"
    printMatrix m4
    
    case invertMatrix m4 of
        Just minv -> do
            putStrLn "\nInverse matrix:"
            printMatrix minv
            
        Nothing -> putStrLn "Matrix is not invertible!"
