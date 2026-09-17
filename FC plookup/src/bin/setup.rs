use fc_plookup::MAX_DEGREE;
use plookup::kzg10;

/// KZG10 的 trusted setup 內部用固定種子的 `test_rng()`，
/// 所以 setup / prover / verifier 各自呼叫 `kzg10::trusted_setup` 就能得到完全相同的參數，
/// 不需要（也沒辦法簡單地）把 `Powers` 序列化成檔案再傳遞。
/// 這支程式只是驗證一次可信設定確實能成功產生，並印出參數量級。
fn main() {
    let (powers, _vk) = kzg10::trusted_setup(MAX_DEGREE);

    println!("KZG10 trusted setup 完成（toy protocol，非正式 MPC ceremony）");
    println!("最大多項式次數：{MAX_DEGREE}");
    println!("powers_of_g 數量：{}", powers.powers_of_g.len());
}
