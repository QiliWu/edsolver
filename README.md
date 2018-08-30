# edsolver

EDSolver是一款专注于解决透射电镜数据分析中的电子衍射标定的产品，其目的是让电子衍射分析变得简单，方便，快捷，让任何有需要的分析者都可以快速准确地得到结果。

pdf_crawl.py：从csv文件中提取出物相晶胞参数信息，并保存到json文件（用于后期导入到远程mongodb数据库）。共提取出13万条晶胞参数信息。

ED_caculater.py: 具体标定算法的代码实现。

具体使用请前往网站：http://www.edsolver.com/
![EDSolver](https://github.com/QiliWu/edsolver/blob/master/edsolver.png)
