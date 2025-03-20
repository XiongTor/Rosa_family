## 蔷薇科系统发育相关研究

### 本项目主要包含的内容为研究中所涉及的代码脚本，主要包含以下几点：

1. 数据获取

- 基因组数据抓取
- 地理分布数据抓取
- 性状分布数据抓取
- 染色体数目信息抓取

2. 数据过滤

- 基因组数据过滤
- 地理分布数据过滤

3. 数据处理

- 系统发育分析建树流程
- 多样化分析流程
- 祖先重建分析流程
- 网状进化及 ILS 水平估计分析流程

注：以上内容为项目主要内容，具体内容请参考项目文档。

---

### 以下为本研究的流程图-2025.02

![蔷薇科果实系统发育研究](https://github.com/XiongTor/Rosa_family/blob/main/imag/readme_flowwork/%E8%94%B7%E8%96%87%E7%A7%91%E6%9E%9C%E5%AE%9E%E7%B1%BB%E5%9E%8B%E7%B3%BB%E7%BB%9F%E5%8F%91%E8%82%B2%E7%A0%94%E7%A9%B6.svg)
![蔷薇科多倍化及网状进化研究](https://github.com/XiongTor/Rosa_family/blob/main/imag/readme_flowwork/%E8%94%B7%E8%96%87%E7%A7%91%E7%BD%91%E7%8A%B6%E8%BF%9B%E5%8C%96%E5%8F%8A%E4%B8%8D%E5%AE%8C%E5%85%A8%E8%B0%B1%E7%B3%BB%E7%AD%9B%E9%80%89%E7%A0%94%E7%A9%B6.svg)
![蔷薇科核型演化研究](https://github.com/XiongTor/Rosa_family/blob/main/imag/readme_flowwork/%E8%94%B7%E8%96%87%E7%A7%91%E6%A0%B8%E5%9E%8B%E6%BC%94%E5%8C%96.svg)

---
## Scripts
>.exe means the script can be run directly in the command line.

- ### **DATA_DOWNLOAD**
  #### **character data download**
     - **get_character_data.R**  
     
       This script is used to download the chatacter data from the TRY, BIEN,GIFT. And then to combine them into one file.

     - **get_dispersal_type.R**

       This script is used to classify the dispersal type of each species based on the dispersal terms provided in the dataset. And to explore the relashionship between dispersal type and distribution range.

  #### **distribution data download and simple plot**
    - **get_distribution_data.R**

      This script is used to download the distribution data from the GBIF, idigbio, NSII, CVH database.

    - **clean_distribution_data.R (.exe)**

      This script is used to clean the distribution data. The script to be cleaned must contain the **Species**, **Longitude**, and **Latitude** fields as column names.   

      ```{eval=F,echo= F,bash}
      # Example: 
      Rscript clean_distribution_data.R data/distribution_data.csv
      ```

    - **change_longtitude_to_address_baidu.py**

      This script is used to convert the longitude to the corresponding address using the Baidu API.

    - **get_wcvp_native_distribution.R**

      This script is used to combine different datasets to be one dataset. And then to change the species name to the accpeted name. Then to get the native distribution by WCVP.

    - **map_rosa_distribution.R**

      This script is used to map the Rosa subfamilies coordinate information to the world map.

    - **count_continent_num_by_coordinate.R**

      This script is used to count the number of continents by the coordinate information.

---
## Scripts
>.exe means the script can be run directly in the command line.