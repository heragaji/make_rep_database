# make_rep_database

program of detecting and clustering repeat sequences from reads.

## Description

`make_rep_database` detects and classifies the repeat sequences from reads in FASTA file, and outputs the repeat databases with FASTA format.

## Requirement

- seqkit Version: 0.8.0
- samtools Version: 1.7 Using htslib 1.7
- minialign Version: 0.6.0-42-g4d6b274, Build: SSE4.1
- snakemake Version: 5.4.0

realignerとdump-consensusの設定は、ラボ内DocBase参照。

- realigner commit-number: d26f40970f705fde4167ef119402df3330e82fb2

- dump-consensus commit-number: a7733b5e0d89c7536f72c3f4e1d5a0c7747d013f

- multialn  comimt-number: 8102e5be6e1319fb07917c44cb484f564ddd91ff

- haskell-sam   commit-number: e8d1bc609eae5415755ca3084b308dcba12d4d7d

- fastlinkedlist    commit-number: bac809b90ed06b2fbb28b78da7c69021e3641d53

- bioseq    commit-number: d36162340f88169a137797124e437e44a6e688ac

- attoparsec-applicative    comimt-number: bd4b2902570dec6af0078171b94075f25f01e708

- python3

    - pysam　Version: 0.11.2.2

    - wiggelen Version: 0.4.1

    - scipy Version: 1.0.1

    - matplotlib Version: 2.2.2

    - biopython Version: 1.71


## Usage

以下のコマンド。

`snakemake --snakefile [snakefile_path] --directory [working_directory_path]`

作業ディレクトリに、`config.json`という名前のJSONファイルを作成する。  
JSONファイルの中身は、

``` json
{
    "data": "入力のFASTAファイル(絶対パス)",
    "output": "出力ファイルを出すディレクトリを指定(絶対パス)",
    "realigner": "realignerのディレクトリ(絶対パス)",
    "dump_consensus": "dump-consensusのディレクトリ(絶対パス)",
    "top": "長い順に上から何本をrepeat検出のreferenceとして使うか",
    "iteration": "realignerを最大何回回すか",
    "coverage": "カバレッジいくつ以上の領域をリピートとみなすか",
    "interval": "何bp内のリードの終端のカバレッジを極大値とみなすか",
    "peak": "リードの終端の数がいくつ以上のとき、ピークとみなすか",
    "cut": "最小リピートの長さ(bp)",
}
```

出力として、"output"で指定したディレクトリに、topで使われた各リードの名前のフォルダ(中間ファイルが入っている), top.fa, repeat.fa,top_vs_reads_sorted.bam,rep_vs_rep.paf,strong.paf,weak.paf,read_masked.faが入る。この中のrepeat.faがrepeatのデータベースとなる。
作業ディレクトリには.snakemakeフォルダ(snakemakeのログ等が入る)が入る。
`--cores int`,`--qsub int`等でコア数とジョブの投入のオプションもできる。

## Installation

    $ git clone https://github.com/heragaji/make_rep_database.git

## Author

Taro Ozawa, Yuta Ito

## License
