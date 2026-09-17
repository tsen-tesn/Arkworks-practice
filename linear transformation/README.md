# Arkworks 練習 - linear Transformation

神經網路最基本的單位就是權重（Weight）、輸入（Input）與偏差（Bias），即
$$y = w \cdot x + b$$


## 重點
1.  變數宣告： 如何在 arkworks 中區分 Public Input 與 Private Witness。 (在這個算式中 public input 是 $y$, private witness 是 $w$, $b$。 提問: $x$ 在什麼情況下是 private witness，什麼情況下是 public input? 答案在最後面)
3.  基本運算： 使用 Constraint System 的加法與乘法。
4.  約束檢查： 寫出 enforce_constraint 來確保 $w \cdot x$ 的結果加上 $b$ 真的等於 $y$。

## 提醒
先不做浮點數運算，我們看過 paper 都知道 ZKP 是 on Finite Field 沒辦法使用 IEEE 754 表達 float。如果要做 float 應該要先做 quantiztion (量化) 或是定點數化，所以我們之後再做。
真實世界我們要做 implementation 幾乎都是作浮點數運算，但這是練習所以先嘗試就好

請注意: 寫 ZKP 程式不是像 python 一樣做運算，而是把變數配置到電路中，然後**建立它們之間的等式關係**。

## implementation requirement
$$y = w \cdot x + b$$


請實現以下兩種 y,w,x,b 的組合 (y 請自己算)
1. $w = 3, x = 6, b = 3$
2. $$W = \begin{pmatrix} 2 & 1 \\ 3 & 4 \end{pmatrix}, \quad X = \begin{pmatrix} 5 \\ 2 \end{pmatrix},\quad B = \begin{pmatrix} 3 \\ 1 \end{pmatrix}$$
    提示
    - 請用迴圈動態產生電路約束
    - 需要使用中繼變數 (在 R1CS 中，限制條件必須嚴格遵守 $A \times B = C$。
當我們要計算 $Y_0 = W_{00} \cdot X_0 + W_{01} \cdot X_1$ 時，因為 $W$ 和 $X$ 都是未知數（Witness），沒辦法把兩個乘法項加在一起直接塞進一個約束裡。需要宣告中繼變數（Intermediate Witness）來儲存乘法結果。)
    
    
3. 請完成下列任務
    - 核心電路實作 (```src/lib.rs```)
        - 建立兩個的結構體 (就是面寫的兩個$ywxb$組合)，並實作結構體的ConstraintSynthesizer。
        - 提示: ```LinearCombination``` 處理加法(加法cost免費)，使用中繼變數。
    - 可信設定 (```src/bin/setup.rs```)
        - 使用 ```ark_bls12_381::Fr``` 作為 Finite Field
        - 使用  ```ark_groth16::Groth16 ``` 搭配一個隨機數生成器（例如  ```test_rng() ```）進行  ```generate_random_parameters_with_reduction ```。
        - 將生成的 pk 與 vk 序列化，並分別匯出成實體檔案： ```pk.bin``` 與  ```vk.bin```。
    - 證明生成與儲存 (```src/bin/prover.rs```)
        - 讀取 ```pk.bin```
        - 使用```Groth16::prove```生成證明
        - 序列化要求：請將產生的 Proof 分別用兩種方式存檔:
            - proof_compressed.data（使用 serialize_compressed）
            - proof_uncompressed.data（使用 serialize_uncompressed）
            - 兩者差異是**大小**，下面那個是 verifier 不需在額外算東西，上面要
            
    - 獨立驗證程式 (```src/bin/verifier.rs```)
        - 這是一個完全獨立的執行檔。不要在這個檔案中實作 LinearLayerCircuit 結構。
        - 從磁碟讀取 vk.bin，並呼叫 ark_groth16::prepare_verifying_key 進行預先計算。
        - 從磁碟讀取 proof_compressed.data 或 proof_uncompressed.data。
        - 在程式碼中手動宣告 Public Input。
        - 呼叫 Groth16::verify_with_processed_vk 進行驗證，並在終端機印出驗證結果 (True/False)。
        - (opt.) 除錯測試：嘗試在驗證時，傳入錯誤的 Public Input（例如 $Y = [15, 25]$），確認 Verifier 能夠回傳 False。



4. terminal 指令
    - ```cargo run --bin setup```
    - ```cargo run --bin prover```
    - ```cargo run --bin verifier```




### Answer of $x$
$x$ 如果以 server-client 角度來說 (要保護 server 隱私的話)，$x$ 其實應該是 user 的 input 所以應該是公開輸入。

還有另一種情境: 要保護client 隱私的話
假設 Server 提供了一個公開的模型（$W, B$），或者是把權重下發給了 Client。Client 想用自己的敏感資料（例如個人的醫療數據、金融資產 $x$）跑模型得到結果 $y$，並向 Server 證明「我算出來的結果符合某種資格」，但絕對不能把真正的 $x$ 傳給 Server。 (by AI)
這樣的情況 $x$ 就是私有輸入

在這個實作中你們可以隨便選擇 應該是讓$x$是 private 會比較好做，請自行選擇 或是都嘗試看看(可能很花時間，不強求)。