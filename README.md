# Oephys_analysis-tools
Neuropixels spinal cord recording using Openephys

## 使用方法
1. 全てのファイルにパスを通す<br>
2. analysis-tools_originalフォルダ内のopen_ephys_disp_exp.mを実行する<br>
3. GUIが表示されるので、"Setting"の項目で、以下の条件を指定する<br>
 ・PC：解析に使用しているPCのOS<br>
 ・Probe：複数本使用している場合の、解析対象プローブ（1本ずつしか解析できません）<br>
 ・Number of culumns: Neuropixelsの列数（4列1bankか、2列2bankか）<br>
4. Purposeの項目で、"LFP mapping"を選択する<br>
5. OKをクリックする<br>
6. 別のGUIが表示されるので、以下の解析のパラメータを入力する<br>
 ・Trigger channel: 刺激装置のトリガー信号が入力されているAIの番号<br>
 ・CDP channel: 銀ボール電極で記録している場合、その信号が入力されているAIの番号<br>
 ・Spike detection channel: スパイク検出結果を表示するチャンネル番号<br>
 ・Nuber of triggers：加算平均回数　多いほど解析に時間を要する<br>
 ・Gain: 解析結果を表示する際の縦軸　大きくすると波形が大きく表示される　大きすぎるようなら小さくすると良い<br>
7. OKをクリックする<br>
8. 解析結果が表示される　結果のFigureは一定時間が経過すると画面からは消えるが、自動的に保存されている<br>
