# FC_plookup 實作設計文件

> 本檔案補充 `README.md` 沒有講清楚的部分：**如何在保持 $z$ 私密的前提下，把「$z=WX+b$ 的算術正確性」跟「$(z_i,y_i)$ 落在 ReLU table 裡」這兩段證明綁在一起**。README 只畫出資料流（$W,X,b \to z \to y$），沒有處理「兩段證明各自驗證，沒有東西保證雙方用的是同一個 $z$」這個問題。這份文件是可以照著實作的步驟拆解。

---

## 0. 問題重述：為什麼「兩段分開證明」不夠

把整個系統拆成：

- **證明 A**（算術）：「我知道 $W,X,b$，使得 $z=WX+b$」
- **證明 B**（plookup）：「我手上的 $(z_i,y_i)$ 都在 ReLU table 裡」

各自獨立驗證的話，惡意 prover 可以：

1. 用真實 $W,X,b$ 算出真正的 $z=[5,-3,2]$，產生合法的證明 A。
2. 在證明 B 完全換一組跟 A 無關、但一樣查表會過的假資料，例如 $z'=[100,100,100],\,y'=[100,100,100]$（因為 $\mathrm{ReLU}(100)=100$）。

兩個證明各自驗證都會通過，但「$y=\mathrm{ReLU}(WX+b)$」這句話其實是假的——這就是 soundness 被破壞的具體攻擊。

**這份文件要解決的唯一問題：讓 A 跟 B 密碼學上被強迫使用同一個 $z$，而且過程中 $z$ 不曾被揭露給 verifier。**

### 為什麼不能簡單把 $z$ 設成 public

如果讓 $z$ 公開，兩段證明自然就綁定了（verifier 直接拿同一組數字餵給兩邊），但這樣做等於放棄了使用 plookup 的理由：

- Plookup 存在的意義就是「證明一個**私密** witness 落在公開表裡，且不用揭露它」。
- 若 $z$ 已公開，verifier 自己就能算 $\mathrm{ReLU}(z)$ 跟 $y$ 比對，完全不需要密碼學查表機制。

所以正確做法必須讓 $z$ 全程保密，同時仍然被綁定。

### 為什麼不能只靠 KZG commitment 本身

KZG commitment 只有**加法同態性**：

$$\mathrm{commit}(a\cdot p(X)+b\cdot q(X)) = a\cdot\mathrm{commit}(p)+b\cdot\mathrm{commit}(q)\quad(a,b\text{ 為公開純量})$$

但 $w_{ij}\cdot x_j$ 是兩個**都保密**的值相乘，KZG 沒辦法只靠 commitment 驗證這種乘法關係——這正是當初要發明 R1CS / PLONK 這類電路證明系統的根本原因。**沒有任何方式可以跳過電路化，直接靠幾個 commitment 檢查乘積。**

### 為什麼選擇「手刻 mini-PLONK」而不是接 Groth16

`linear transformation/` 現有的矩陣乘法是用 `ark_relations` 的 R1CS + Groth16，這是一套跟 `plookup/` 完全不同的密碼學世界（pairing-based SNARK vs. 自製的 KZG10 + multiset argument）。硬接兩者需要額外一層「跨系統」的綁定機制，而且 Groth16 的 proof 結構沒辦法讓 verifier 用同態性質自己組出 commitment。

`plookup/` 現有的程式碼（`kzg10.rs`、`transcript.rs`、`quotient_poly.rs`、`multiset_equality.rs`）其實已經是一個小型 PLONK 工具箱：FFT domain、KZG commit/open、Fiat-Shamir transcript、grand-product accumulator、quotient polynomial 除以 vanishing polynomial 這些機制全部都在。**只要把 $W,X,b$ 的乘法/加法關係也用同一套語言（gate + selector + wire polynomial）表達，算術證明跟 lookup 證明就會天生活在同一個框架裡**，不需要另外發明綁定機制。

---

## 1. 整體架構總覽

```text
                 ┌─────────────────────────────────────────┐
                 │   Mini-PLONK 電路（M1~M3）                │
                 │   wire: a(X), b(X), c(X)                 │
                 │   selector: q_M,q_L,q_R,q_O,q_C（公開）   │
                 │                                          │
  W, X, b ─────▶│  M2: Gate identity（乘法/加法正確性）      │
   （私密）       │  M3: Permutation argument（wiring 正確性）│
                 │                                          │
                 │  輸出：z1,z2,z3 是 c(X) 在固定位置的值      │
                 └───────────────┬──────────────────────────┘
                                 │ C_z（KZG commitment，z 全程不揭露）
                                 ▼
                 ┌─────────────────────────────────────────┐
                 │   既有 plookup 模組（略修改）                │
                 │   M4: f_commit = C_z + θ·C_y（同態組出）    │
                 │   （沿用 multiset_equality + quotient_poly）│
                 └───────────────┬───────────────────────────┘
                                 │
                                 ▼
                 M5: 跨 domain 開啟一致性證明（真正鎖住 z 的地方）
                                 │
                                 ▼
                 M6: 全程共用一條 Fiat-Shamir transcript
```

---

## 2. M1：Wire 佈局與 Gate 表（3×3 矩陣乘法攤平）

把 $z_i=\sum_j w_{ij}x_j+b_i$ 攤成基本閘。每列 4 項相加要拆成 3 個二元加法：

$$t_{i1}=p_{i0}+p_{i1},\quad t_{i2}=t_{i1}+p_{i2},\quad z_i=t_{i2}+b_i$$

共 $9$ 個乘法閘 $+9$ 個加法閘 $=18$ 個閘。每個閘是一個 $(a_k,b_k,c_k)$ 三元組，補到 $N=32$（$2$ 的冪次，多出的閘填 dummy `0+0=0`，用跟真實乘法閘相同的 selector，只是輸入輸出都填 0，恆等式自動成立）。

### 2.1 乘法閘（gate index $0\sim8$）

| gate $k$ | $a_k$ | $b_k$ | $c_k$ | selector |
|---:|---|---|---|---|
| 0 | $w_{00}$ | $x_0$ | $p_{00}$ | $q_M=1,q_O=-1$ |
| 1 | $w_{01}$ | $x_1$ | $p_{01}$ | 同上 |
| 2 | $w_{02}$ | $x_2$ | $p_{02}$ | 同上 |
| 3 | $w_{10}$ | $x_0$ | $p_{10}$ | 同上 |
| 4 | $w_{11}$ | $x_1$ | $p_{11}$ | 同上 |
| 5 | $w_{12}$ | $x_2$ | $p_{12}$ | 同上 |
| 6 | $w_{20}$ | $x_0$ | $p_{20}$ | 同上 |
| 7 | $w_{21}$ | $x_1$ | $p_{21}$ | 同上 |
| 8 | $w_{22}$ | $x_2$ | $p_{22}$ | 同上 |

恆等式：$a_k b_k - c_k = 0$。

### 2.2 加法閘（gate index $9\sim17$）

| gate $k$ | $a_k$ | $b_k$ | $c_k$ | selector |
|---:|---|---|---|---|
| 9  | $p_{00}$（=$c_0$） | $p_{01}$（=$c_1$） | $t_{01}$ | $q_L=1,q_R=1,q_O=-1$ |
| 10 | $t_{01}$（=$c_9$） | $p_{02}$（=$c_2$） | $t_{02}$ | 同上 |
| 11 | $t_{02}$（=$c_{10}$） | $b_0$ | $z_0$ | 同上 |
| 12 | $p_{10}$（=$c_3$） | $p_{11}$（=$c_4$） | $t_{11}$ | 同上 |
| 13 | $t_{11}$（=$c_{12}$） | $p_{12}$（=$c_5$） | $t_{12}$ | 同上 |
| 14 | $t_{12}$（=$c_{13}$） | $b_1$ | $z_1$ | 同上 |
| 15 | $p_{20}$（=$c_6$） | $p_{21}$（=$c_7$） | $t_{21}$ | 同上 |
| 16 | $t_{21}$（=$c_{15}$） | $p_{22}$（=$c_8$） | $t_{22}$ | 同上 |
| 17 | $t_{22}$（=$c_{16}$） | $b_2$ | $z_2$ | 同上 |

恆等式：$a_k+b_k-c_k=0$。

**$z_1,z_2,z_3$ 分別是 gate 11、14、17 的輸出 $c$ 值**，也就是 $c(\omega^{11}),c(\omega^{14}),c(\omega^{17})$。這三個 domain 點的位置是公開、固定的（電路拓樸本身不是秘密，只有填進去的值是秘密）。

---

## 3. M2：Gate Identity（算術正確性）

把 selector 向量跟 wire 值向量分別透過 `domain.ifft`（跟 `MultiSet::to_polynomial` 現在做的事一模一樣）轉成多項式 $q_M(X),q_L(X),q_R(X),q_O(X),q_C(X),a(X),b(X),c(X)$。

構造：

$$\mathrm{GateCheck}(X)=q_M(X)a(X)b(X)+q_L(X)a(X)+q_R(X)b(X)+q_O(X)c(X)+q_C(X)$$

若每個閘都被正確滿足，這個多項式在整個 $32$ 點 domain 上每一點都是 $0$，也就是能被 vanishing polynomial $Z_H(X)=X^{32}-1$ 整除，餘式為零。

實作上完全比照 `plookup/src/multiset/quotient_poly.rs::compute` 現在的模式：

```rust
let gate_check = &(&(&q_m * &(&a * &b)) + &(&q_l * &a)) + &(&(&q_r * &b) + &(&q_o * &c)) + &q_c;
let (quotient, remainder) = gate_check.divide_by_vanishing_poly(domain);
assert!(remainder.is_zero()); // 合法 witness 時必須成立
```

> **建議先做的最小可測試單元**：只做 M1+M2，寫一個獨立函式，輸入合法/非法的 $w,x,b$，輸出 remainder 是否為零。這一步不依賴 M3～M6，可以先驗證電路佈局思路正確。

---

## 4. M3：Permutation Argument（wiring 正確性）

Gate identity 只保證「每個閘自己內部算對了」，不保證「乘法閘的輸出真的被接到對的加法閘輸入」。這需要 copy constraint（PLONK 的 permutation argument）。

### 4.1 需要哪些 copy constraint

**乘法 → 加法的接線**（每組是一個 2-cycle）：

| 來源 cell | 目標 cell |
|---|---|
| $c_0$（$p_{00}$） | $a_9$ |
| $c_1$（$p_{01}$） | $b_9$ |
| $c_9$（$t_{01}$） | $a_{10}$ |
| $c_2$（$p_{02}$） | $b_{10}$ |
| $c_{10}$（$t_{02}$） | $a_{11}$ |
| $c_3$（$p_{10}$） | $a_{12}$ |
| $c_4$（$p_{11}$） | $b_{12}$ |
| $c_{12}$（$t_{11}$） | $a_{13}$ |
| $c_5$（$p_{12}$） | $b_{13}$ |
| $c_{13}$（$t_{12}$） | $a_{14}$ |
| $c_6$（$p_{20}$） | $a_{15}$ |
| $c_7$（$p_{21}$） | $b_{15}$ |
| $c_{15}$（$t_{21}$） | $a_{16}$ |
| $c_8$（$p_{22}$） | $b_{16}$ |
| $c_{16}$（$t_{22}$） | $a_{17}$ |

**$x_j$ 被重複使用**（每個是一個 3-cycle，因為同一個 $x_j$ 在 3 個列的乘法閘裡都當 $b$ 輸入）：

- $\{b_0, b_3, b_6\}$ 都必須等於 $x_0$
- $\{b_1, b_4, b_7\}$ 都必須等於 $x_1$
- $\{b_2, b_5, b_8\}$ 都必須等於 $x_2$

$w_{ij}$ 跟 $b_i$（bias）各自只出現一次，不需要額外 copy constraint（在置換 $\sigma$ 裡是自身對應的 fixed point）。

### 4.2 Grand-product accumulator

這一步跟 `multiset_equality.rs::compute_accumulator_values` **完全同構**：一樣是 grand-product accumulator $Z(X)$、一樣有 $L_0,L_{n}$ 邊界條件、一樣有一條 transition identity 要除以 vanishing polynomial。差別只在乘積項換成標準 PLONK permutation check：

$$Z(\omega X)\cdot\prod_{\text{wire columns}}\big(\text{cell}+\beta\cdot\sigma(\text{id})+\gamma\big) \;=\; Z(X)\cdot\prod_{\text{wire columns}}\big(\text{cell}+\beta\cdot\text{id}+\gamma\big)$$

其中 `id` 是每個 cell 的原始編號（依 column 跟 domain 位置），`σ(id)` 是上面表格定義出的置換結果，$\beta,\gamma$ 是 Fiat-Shamir challenge。

> 這是整個 mini-PLONK 裡最繁瑣的部分。建議直接對照 PLONK 論文第 5 節的 permutation argument，或參考 `multiset_equality.rs` + `quotient_poly.rs` 的 `compute_term_check_a/b` 寫法去改寫成 permutation 版本。

---

## 5. M4：把 $z,y$ 接進既有的 plookup 模組——用 KZG 加法同態

目前 `plookup/src/multiset/proof.rs::EqualityProof::prove` 是自己重新 commit 一份 `f_poly`：

```rust
let f_commit = kzg10::commit(proving_key, &f_poly);
```

**改成不要讓 prover 自己重新 commit f**，而是讓 verifier 自己用同態性質算出來：

$$C_f = C_z + \theta\cdot C_y$$

（$C_z$ 是 M1～M3 電路裡 $z$-wire 的 commitment，$C_y$ 是 $y$-wire 的 commitment，$\theta$ 是 Fiat-Shamir challenge，對應現在 `compress_column` 用的同一個 $\theta$）。

因為 KZG commitment 對線性組合是同態的，verifier 完全不需要相信 prover 講的「這是 f 的 commitment」——這份 commitment 是 verifier 自己從已經在 M1～M3 驗證過的 $C_z,C_y$ 算出來的，天生保證 f 用的就是同一個 $z$。

---

## 6. M5：跨 Domain 的開啟一致性（真正鎖住 z 的地方）

問題：M1～M3 的 gate 電路 domain 大小是 $N=32$，但既有 `plookup` 模組內部工作 domain 綁在 `TABLE_SIZE=64`。$C_z$（32-domain 上 $z$-wire 的 commitment）沒辦法直接套用在 M4 的線性組合，因為那是不同 domain 上的多項式。

修法：額外加一組（可批次聚合成一次）KZG 單點開啟證明：

1. 對 $C_z$（32-domain 的 wire commitment）在固定點 $\omega^{11},\omega^{14},\omega^{17}$ 各開一個 KZG 單點證明，開出 $z_0,z_1,z_2$。
2. 對 plookup 內部的 $C_f$（64-domain）在 $z_i$ 對應的 witness 位置也開一個單點證明，並用 `evaluation = z_i + \theta\cdot y_i` 反推出同一個 $z_i$。
3. Verifier 檢查兩邊開出來的 $z_i$ 相等。

三個點可以比照 `plookup/src/multiset/proof.rs` 裡 `compute_aggregate_witness` 的做法，用一個 aggregation challenge 打包成一次批次開啟證明，不用開 3 次。

**這一步就是整個設計裡「真正讓 A、B 兩段證明用同一個 z」的地方**：M4 保證 f 是 $z,y$ 的線性組合，M5 保證這裡的 $z$ 跟電路算出來、經過 gate identity + permutation 驗證過的 $z$ 是同一個數字，而且全程沒有任何一步把 $z$ 的實際數值透露給 verifier。

---

## 7. M6：只能有一條 Fiat-Shamir Transcript

M1～M5 的所有 commitment（gate wire、permutation 用的 $\beta,\gamma$、plookup 的 $h_1,h_2,Z,$ quotient commitment、M5 的開啟證明挑戰）都必須**依序**餵進同一個 `merlin::Transcript`，順序寫死、不可調換。這樣才是一個不可分割的 non-interactive proof，而不是兩個各自做 Fiat-Shamir 的獨立 proof 硬貼在一起（否則 binding 問題只是換了個位置重新出現：例如 permutation 的 $\beta,\gamma$ 若用獨立 transcript 產生，理論上仍可能被抓到操作空間）。

建議的 transcript label 順序：

```text
gate_a_commit, gate_b_commit, gate_c_commit
  → beta_perm, gamma_perm（M3 challenge）
perm_z_commit（M3 accumulator）
gate_quotient_commit（M2+M3 合併的 quotient）
  → evaluation_challenge（M2/M3 的開值點）
z_wire_eval, y_wire_eval, ...（M5 開值）
theta（M4 的 lookup 壓縮 challenge，注意要晚於 z,y 的 commitment 才產生，避免舊有「theta 寫死成常數」的問題一併修正）
h_1_commit, h_2_commit
  → beta_lookup, gamma_lookup
accumulator_commit（plookup 的 Z(X)）
quotient_commit（plookup 的 quotient）
  → evaluation_challenge_lookup
（各 evaluation 值）
  → aggregation_challenge
```

---

## 8. 建議實作順序 / 檔案結構

1. `src/gate/layout.rs`：寫死 2.1、2.2 節的 gate 表（純資料，no crypto），提供「輸入 $w,x,b$ → 輸出 $a,b,c$ 三條 wire 的完整向量」的函式。先寫單元測試確認佈局本身算出來的 $z_i$ 是對的。
2. `src/gate/quotient.rs`：實作 M2 的 `GateCheck` + `divide_by_vanishing_poly`，先只測「remainder 是否為零」，不涉及 commitment。
3. `src/gate/permutation.rs`：實作 M3，抄 `multiset_equality.rs` 的 accumulator 結構，換成 permutation 版本的乘積項。
4. `src/gate/proof.rs`：把 M2+M3 包成一個完整的 prove/verify（commit wire、commit selector 可省略因為公開、commit accumulator、commit quotient、開值、batch verify），模仿 `plookup/src/multiset/proof.rs::EqualityProof` 的結構。
5. 修改（或包一層）`plookup/src/multiset/proof.rs`：加入 M4 的同態組合 `f_commit`、拿掉現在寫死的 `theta = Fr::from(7u64)`，改成 transcript challenge。
6. `src/binding.rs`：實作 M5 的跨 domain 開啟一致性證明。
7. `src/fc_proof.rs`：串起 M1～M6，提供最外層 `prove(w,x,b) -> FCProof` / `verify(proof, table_commitment) -> bool`，內部維護單一 `Transcript`。

每一步都建議先寫小型單元測試（用具體數字，例如 README 範例的 $W,X,b$）確認该步驟本身的 remainder / accumulator 邊界條件正確，再往下一步疊加，避免最後才整合、debug 範圍過大。
