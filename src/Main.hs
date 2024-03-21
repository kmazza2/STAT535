module Main where
sqPlusTwo x = x * x + 2
main :: IO ()
main = do
          let z = sqPlusTwo 3
          putStrLn $ "The result is: " ++ show z
