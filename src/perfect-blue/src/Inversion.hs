{-#Language UndecidableInstances,ScopedTypeVariables,RankNTypes,TypeOperators,DataKinds,StandaloneDeriving,KindSignatures,TypeApplications,MultiParamTypeClasses 
           ,GADTs,FlexibleInstances,ConstraintKinds,FunctionalDependencies,FlexibleContexts,AllowAmbiguousTypes#-}

module Inversion where
-- Pure Haskell matrix inversion using Gauss-Jordan elimination
-- No external BLAS/LAPACK dependencies

type Matrix = [[Double]]

-- Clean up near-zero values
cleanZeros :: Matrix -> Matrix
cleanZeros = map (map clean)
  where clean x = if abs x < 1e-10 then 0 else x

-- Pretty print a matrix
printMatrix :: Matrix -> IO ()
printMatrix m = mapM_ printRow m
  where
    printRow row = putStrLn $ unwords $ map (printf "%.4f") row
    printf fmt x = take 10 $ show x ++ replicate 10 ' '

-- Matrix inversion using Gauss-Jordan elimination
-- Returns Nothing if the matrix is singular
invertMatrix :: Matrix -> Maybe Matrix
invertMatrix m
    | rows /= cols = Nothing
    | otherwise = gaussJordan augmented
  where
    rows = length m
    cols = length (head m)
    identity = [[ if i == j then 1 else 0 | j <- [0..cols-1]] | i <- [0..rows-1]]
    augmented = zipWith (++) m identity

-- Gauss-Jordan elimination
gaussJordan :: Matrix -> Maybe Matrix
gaussJordan m = do
    reduced <- forwardElimination m 0
    result <- backSubstitution reduced (length reduced - 1)
    return $ map (drop (length result)) result

-- Forward elimination to get row echelon form
forwardElimination :: Matrix -> Int -> Maybe Matrix
forwardElimination m row
    | row >= length m = Just m
    | otherwise = do
        m' <- pivotRow m row
        let m'' = eliminateColumn m' row
        forwardElimination m'' (row + 1)

-- Find and swap pivot row
pivotRow :: Matrix -> Int -> Maybe Matrix
pivotRow m row =
    case findPivot m row row of
        Nothing -> Nothing
        Just pivotIdx -> Just $ swapRows m row pivotIdx

-- Find non-zero pivot
findPivot :: Matrix -> Int -> Int -> Maybe Int
findPivot m col row
    | row >= length m = Nothing
    | abs (m !! row !! col) > 1e-10 = Just row
    | otherwise = findPivot m col (row + 1)

-- Swap two rows
swapRows :: Matrix -> Int -> Int -> Matrix
swapRows m i j =
    [ if k == i then m !! j
      else if k == j then m !! i
      else m !! k
    | k <- [0..length m - 1]]

-- Eliminate column below pivot
eliminateColumn :: Matrix -> Int -> Matrix
eliminateColumn m row =
    [ if i <= row then m !! i
      else let pivot = m !! row !! row
               factor = (m !! i !! row) / pivot
           in zipWith (\a b -> a - factor * b) (m !! i) (m !! row)
    | i <- [0..length m - 1]]

-- Back substitution to get reduced row echelon form
backSubstitution :: Matrix -> Int -> Maybe Matrix
backSubstitution m row
    | row < 0 = Just m
    | abs (m !! row !! row) < 1e-10 = Nothing
    | otherwise = do
        let m' = scaleRow m row
            m'' = eliminateAbove m' row
        backSubstitution m'' (row - 1)

-- Scale row so diagonal element is 1
scaleRow :: Matrix -> Int -> Matrix
scaleRow m row =
    [ if i /= row then m !! i
      else let pivot = m !! row !! row
           in map (/ pivot) (m !! i)
    | i <- [0..length m - 1]]

-- Eliminate column above pivot
eliminateAbove :: Matrix -> Int -> Matrix
eliminateAbove m row =
    [ if i >= row then m !! i
      else let factor = m !! i !! row
           in zipWith (\a b -> a - factor * b) (m !! i) (m !! row)
    | i <- [0..length m - 1]]

-- Matrix multiplication
matMul :: Matrix -> Matrix -> Matrix
matMul a b =
    [[sum $ zipWith (*) row col | col <- transpose b] | row <- a]
  where
    transpose = foldr (zipWith (:)) (repeat [])

main :: IO ()
main = do
    -- Create a 3x3 invertible matrix
    let m = [[4, 7, 2],
             [1, 3, 1],
             [2, 5, 3]] :: Matrix
    
    putStrLn "Original matrix:"
    printMatrix m
    
    case invertMatrix m of
        Just minv -> do
            putStrLn "\nInverse matrix:"
            printMatrix minv
            
            -- Verify: M * M^-1 should be identity
            let identity = cleanZeros $ matMul m minv
            putStrLn "\nVerification (M * M^-1):"
            printMatrix identity
            
        Nothing -> putStrLn "Matrix is not invertible!"
    
    -- Test with singular matrix
    putStrLn "\n--- Testing singular matrix ---"
    let singular = [[1, 2, 3],
                    [2, 4, 6],
                    [1, 2, 3]] :: Matrix
    
    putStrLn "Singular matrix:"
    printMatrix singular
    
    case invertMatrix singular of
        Just _ -> putStrLn "Unexpectedly got an inverse"
        Nothing -> putStrLn "\nMatrix is not invertible (as expected)"