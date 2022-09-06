# Astro Image Measurement Tool v0.1

## 概要

`astro-image-measurement.py` は、天体写真に写った星の画像上の座標からその星の位置(赤経/赤緯)を計算してアノテーション付きの画像を生成するツールです。

## 注意

このツールは天文学の教育を受けたことのないアマチュア天文家が個人的な趣味で作成したもので、計算結果の正確性は保証できません。天文学の専門家のレビューを受けるまで 1.0 未満のバージョンでリリースする予定です。

## 動作環境

Python 3 が動く環境で動作します。
以下の環境で動作確認しています。

- Ubuntu 20.04 LTS
- Anaconda 4.12
  - python 3.10.4
  - astropy 5.0.4
  - astroquery 0.4.6
  - svgwrite 1.4.2
  - pillow 9.0.1

## サンプル

入力/出力のサンプルファイルが [sample](/sample) フォルダに入っています。

- [in_image.png](/sample/in_image.png): 入力画像(シリウス周辺を撮影したもので中央部は短時間露出で撮ったものを合成しています。)
- [reference_stars.json](/sample/reference_stars.json): 位置決めの基準となる星の名前と座標を記述した入力ファイル。
- [target_stars.json](/sample/target_stars.json): 位置を知りたい星の画像上での座標、または位置が既知で画像上での座標を知りたい星の名前を記述した入力ファイル。
- [out_image.svg](/sample/out_image.svg): in_image.png の画像上に入力ファイルに記述された星のアノテーションと赤経/赤緯のグリッド線を描画した出力ファイル。
- [result.json](/sample/result.json): 計算結果が記述された出力ファイル。
- [out_image.png](/sample/out_image.png): out_image.svg を Incscape で PNG 画像に変換したものです。SVG 画像の描画の確認用です。出力画像の SVG ファィルの編集には、この画像と同じように表示できるツールを使用してください(フォントの表示は環境によって異なります)。

## 使用法

### 参照星データファイルの作成

参照星データファイル(reference stars data file)は、画像上の座標系と天空座標系の対応関係を決めるための基準となる星(参照星)の名前と座標を記述した入力ファイルです。サンプルでは [reference_stars.json](/sample/reference_stars.json) がそれです。

参照星データファイルは、以下の手順で作成します。

1. 入力画像上から [SIMBAD](http://simbad.cds.unistra.fr/simbad/) に登録のある既知の恒星を3つ以上特定し、画像上の座標を測定します。
   - 画像の左上の隅が原点の座標系で測定します。
   - 単位はピクセル単位です。
2. 参照星データファイルに、以下を記述します。
   - 入力画像の撮影日時。
   - 1で特定した星の名前と座標。

以下が参照星データファイルの例です。

```json
{
  "obstime": "2022-03-12 10:52:07",
  "stars": [
    {
      "name": "BD-16 1589",
      "x": 3826,
      "y": 3398
    },
    {
      "name": "BD-16 1586",
      "x": 5149,
      "y": 2276
    },
    {
      "name": "Gaia DR2 2947057475918406784",
      "x": 2457,
      "y": 888
    }
  ]
}
```

参照星データファイルはJSON形式のファイルで、1つのオブジェクトを記述します。オブジェクトのデータ構造は以下の通りです。

- `obstime`: 入力画像の撮影日時です。世界時(UT)で YYYY-MM-DD hh:mm:ss の書式で記述します。
- `stars`: 参照星のデータを記述したオブジェクトの配列です。
  - 参照星オブジェクト:
    - `name`: 参照星の名前です。SIMBAD に登録された名前(別名も可)を正確に転記します。
    - `x`: 参照星の入力画像上でのx座標値です。単位はピクセルです。
    - `y`: 参照星の入力画像上でのy座標値です。単位はピクセルです。

### 対象星データファイルの作成

対象星データファイル(target stars data file)は、位置を測定したい星の画像上での座標、または既知の星の名前を記述した入力ファイルです。サンプルでは [target_stars.json](/sample/target_stars.json) がそれです。

対象星データファイルは、以下の手順で作成します。

1. 入力画像上から任意の数の位置を測定したい星(対象星)の画像上の座標を測定します。
   - 画像の左上の隅が原点の座標系で測定します。
   - 単位はピクセル単位です。
2. 必要に応じて入力画像上にある参照星以外の任意の数の恒星について SIMBAD に登録された名前を特定します。
3. 対象星データファイルに、以下を記述します。
   - 対象星の名前(仮名でもOK)と1で特定したその座標。
   - 2で特定した星の名前。

以下が対象星データファイルの例です。

```json
[
  {
    "name": "Sirius A"
  },
  {
    "name": "Sirius A (observed (2022/3/12))",
    "x": 3240,
    "y": 2161
  },
  {
    "name": "Sirius B (observed (2022/3/12))",
    "x": 3193,
    "y": 2135
  }
]
```

対象星データファイルはJSON形式のファイルで、任意の数の対象星オブジェクトを含む1つの配列を記述します。対象星オブジェクトのデータ構造は以下の通りです。

- 対象星オブジェクト:
  - `name`: 対象星の名前です。`x`, `y` を記述しない場合は SIMBAD に登録された名前(別名も可)を正確に転記します。`x`, `y` を記述する場合は任意の名前を記述できます。
  - `x`: 対象星の入力画像上でのx座標値です。単位はピクセルです。
  - `y`: 対象星の入力画像上でのy座標値です。単位はピクセルです。

`name` に英語以外(日本語等)の名前を記述する場合はエディタ等で文字コードを UTF-8 に設定して保存してください。

### プログラムの実行

`astro-image-measurement.py` は以下のように実行します。

```sh
python astro-image-measurement.py reference_stars.json target_stars.json in_image.png out_image.svg result.json
```

各引数の意味は以下の通りです。

|引数                      | 意味                                     |
|--------------------------|------------------------------------------|
| reference_stars.json     | 参照星データファイルです。               |
| target_stars.json        | 対象星データファイルです。               |
| in_image.png             | 入力画像ファイルです。JPEG (拡張子 .jpg, .jpeg) と PNG (拡張子 .png) に対応しています。|
| out_image.svg            | 出力画像ファイル名です。SVG 形式のファイル名(拡張子 .svg)を指定します。|
| result.json              | 計算結果データファイルのファイル名です。指定しなかった場合は同じ内容がコンソールに出力されます。|

### 出力画像

出力画像は SVG 形式の画像ファイルです。出力画像は入力画像の上に赤経/赤緯のグリッド線と参照星、対象星のマーカーが描画されたものて、SVG エディタで編集可能です([Inkscape](https://inkscape.org/) での動作を確認しています。ソフトによっては正常に表示・編集できない可能性があります)。

描画オブジェクトはグループ化され、以下のように id が設定されています。

- `original_image`: 入力画像です。元の形式のまま埋め込まれています。
- `grid`: 赤経/赤緯のグリッド線です。
  - `dec`: 赤緯のグリッド線です。
    - `dec_d`: 「度」のグリッド線です。
    - `dec_m`: 「分」のグリッド線です。
    - `dec_s`: 「秒」のグリッド線です。
  - `ra`: 赤経のグリッド線です。
    - `ra_h`: 「時」のグリッド線です。
    - `ra_m`: 「分」のグリッド線です。
    - `ra_s`: 「秒」のグリッド線です。
- `ref_star_1`, `ref_star_2`, ... : 参照星のマーカーです。
  - `ref_star_1_label`, ... : マーカーのラベルです。
    - `ref_star_1_label_name`, ... : 参照星の名前のテキストです。
    - `ref_star_1_label_ra`, ... : 参照星の赤経のテキストです。
    - `ref_star_1_label_dec`, ... : 参照星の赤緯のテキストです。
  - `ref_star_1_crosshair`, ... : マーカーの十字線(星の位置の印)です。
- `target_star_1`, `target_star_2`, ... :
  - `target_star_1_label`, ... : マーカーのラベルです。
    - `target_star_1_label_name`, ... : 対象星の名前のテキストです。
    - `target_star_1_label_ra`, ... : 対象星の赤経のテキストです。
    - `target_star_1_label_dec`, ... : 対象星の赤緯のテキストです。
  - `target_star_1_crosshair`, ... : マーカーの十字線(星の位置の印)です。

各オブジェクトのスタイルは固定です。SVG エディタ等でお好みのスタイルを設定してください。また、近接したマーカーは重なって表示されます。これも SVG エディタ等で手動で調整してください。

### 計算結果データファイル

計算結果データファイルには対象星のデータを計算結果で補完したものです。

以下が計算結果データファイルの例です。

```json
[
  {
    "name": "Sirius A",
    "ra": "6h45m08.07364658s",
    "dec": "-16d43m25.16163781s",
    "x": 3225.698517435161,
    "y": 2157.659906351022
  },
  {
    "name": "Sirius A (observed (2022/3/12))",
    "x": 3240,
    "y": 2161,
    "ra": "6h45m07.85740342s",
    "dec": "-16d43m25.83210444s"
  },
  {
    "name": "Sirius B (observed (2022/3/12))",
    "x": 3193,
    "y": 2135,
    "ra": "6h45m08.57200761s",
    "dec": "-16d43m20.37968445s"
  }
]
```

対象星データファイルはJSON形式のファイルで、計算結果オブジェクトの配列の形をしています。計算結果オブジェクトのデータ構造は以下の通りです。

- 計算結果オブジェクト:
  - `name`: 対象星の名前です。対象星データファイルで記述した名前です。
  - `x`: 対象星の入力画像上でのx座標値です。単位はピクセルです。
  - `y`: 対象星の入力画像上でのy座標値です。単位はピクセルです。
  - `ra`: 対象星の赤経です。
  - `dec`: 対象星の赤緯です。

対象星データファイルで `name` と `x`, `y` を記述した対象星に対応する計算結果オブジェクトでは、`ra`, `dec` が計算結果です。

対象星データファイルで `name` だけを記述した対象星に対応する計算結果オブジェクトでは、`ra`, `dec` が SIMBAD から取得した星の位置データ、`x`, `y` が計算結果で、`ra`, `dec` に対応する入力画像上での座標です。

## 計算方法について

入力画像の座標と天球座標系の相互変換には astropy.wcs.utils の fit_wcs_from_points で生成した WCS を使用しています。

参照星と `x`, `y` を指定しない対象星の位置は、astroquery.simbad で SIMBAD から各星の赤経、赤緯、固有運動、年周視差、視線速度のデータを取得し、そこから astropy で J2000.0 時点の SkyCoord を生成し(SIMBAD から取得した赤経、赤緯を J2000.0 時点のものとしています)、apply_space_motion() で観測時点での位置を計算しています。

## おまけ

`separation.py` を使用すると計算結果データファイルから2つの星の離角を計算できます。

例:

```shell-session
$ python separation.py sample/result.json 'Sirius A (observed (2022/3/12))' 'Sirius B (observed (2022/3/12))'
Sirius A (observed (2022/3/12)): 6h45m07.85740342s / -16d43m25.83210444s
Sirius B (observed (2022/3/12)): 6h45m08.57200761s / -16d43m20.37968445s
separation: 0d00m11.623848s (0d00m10.26567687s / 0d00m05.45241999s)

```

## ライセンス

MITライセンスです。[LICENSE](/LICENSE) を参照してください。
