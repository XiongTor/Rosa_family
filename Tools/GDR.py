from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
import time

# 设置无头模式（可选）
chrome_options = Options()
chrome_options.add_argument("--headless")  # 无头浏览器
chrome_options.add_argument("--disable-gpu")

# 设置你的 chromedriver 路径
driver = webdriver.Chrome(service=Service("chromedriver"), options=chrome_options)
driver.get("https://www.rosaceae.org/tools/jbrowse")
time.sleep(3)  # 等待页面加载，必要！

# 获取所有链接中包含“Genome Page”的链接
genome_links = []
links = driver.find_elements(By.TAG_NAME, "a")
for link in links:
    if "Genome Page" in link.text:
        href = link.get_attribute("href")
        genome_links.append((link.text, href))

print(f"找到 {len(genome_links)} 个 Genome 页面")

# 进入每个 Genome 页面，寻找 “Protein sequences” 下载链接
for name, url in genome_links:
    driver.get(url)
    time.sleep(2)
    try:
        protein_link = driver.find_element(By.PARTIAL_LINK_TEXT, "Protein sequences").get_attribute("href")
        print(f"\n{name}: {url}")
        print(f"  → Protein sequences 下载链接：{protein_link}")
    except:
        print(f"\n{name}: {url}")
        print("  ❌ 未找到 Protein sequences 下载链接")

driver.quit()
