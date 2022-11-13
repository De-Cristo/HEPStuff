1. Varial同时需要ROOT和python2，这是上次和致源调试失败的原因，最简单的办法，CMSSW_10_6_X环境
2. 环境设置还是：https://github.com/De-Cristo/HEPStuff/blob/main/simple_varial/set-varial.sh
3. CPU数量调用：https://github.com/De-Cristo/HEPStuff/blob/main/simple_varial/plot_cfg.py#L10 否则80%全部CPU核数
4. 设置不需要删除原有地址并不需要重新读ntuple：https://github.com/De-Cristo/HEPStuff/blob/main/simple_varial/plot_cfg.py#L14
5. 通过parallel实现ntuple->histo后，通过调用‘lazy function’进行画图或者用已有的histogram做任何工作：https://github.com/De-Cristo/HEPStuff/blob/main/simple_varial/run_plot_from_tree.py#L229-L233
  - note: 在lazy function前修改histograms都会有效应用在后续画图等操作中
  - tool可以是None: https://github.com/De-Cristo/HEPStuff/blob/main/simple_varial/run_plot_from_tree.py#L231
6. 现风格的画图函数在'varial.tools.mk_rootfile_plotter()'，估计（**not been tested**）可以通过修改（或者替换）这个函数，更方便地批量调整histogram（所有操作必须是ROOT中histogram允许的）
