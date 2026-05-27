fn main() {
    // LALRPOPをビルドスクリプトで実行して、parser.rsを生成する
    lalrpop::process_root().unwrap();
}
