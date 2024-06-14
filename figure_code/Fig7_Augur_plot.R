library(ggplot2)
library(ggrepel)

#data_1 <- read.csv("G:\\TLSdata\\TLS_total_data\\ICB_varify\\GSE200996_scRNA\\augur\\Augur_NR.csv",row.names = 1)
#data_2 <- read.csv("G:\\TLSdata\\TLS_total_data\\ICB_varify\\GSE200996_scRNA\\augur\\Augur_R.csv",row.names = 1)

data_1 <- read.csv("G:\\TLSdata\\TLS_total_data\\ICB_varify\\GSE200996_scRNA\\augur\\Augur_NR_single.csv",row.names = 1)
data_2 <- read.csv("G:\\TLSdata\\TLS_total_data\\ICB_varify\\GSE200996_scRNA\\augur\\Augur_R_single.csv",row.names = 1)


data_1 <- data_1[1,]
data_2 <- data_2[1,]

rownames(data_1) <- "NR"
rownames(data_2) <- "R"

data_1 <- as.data.frame(t(data_1))
data_2 <- as.data.frame(t(data_2))

data_1$Key = rownames(data_1)
data_2$Key = rownames(data_2)

data <- merge(data_1, data_2, by = "Key")

data$delta_AUC <- data$R - data$NR

# 为每个点基于其rownames分配颜色
colors <- colorRampPalette(c("#FF0000","#FFFF00","#00FF00",'#00FFFF',"#0000FF","#FF00FF"))(45) ###color
rownames_colors <- setNames(colors, data$Key)


# 计算每个点到直线y=x的垂直距离
distances <- abs(data$R - data$NR) 

# 对距离小于0.1的点将颜色设置为灰色
colors[distances < 0.015] <- "lightgray"



# 绘制散点图并为每个点使用分配的颜色
plot(data$NR, data$R, main = "Augur",
     xlab = "NR", ylab = "R", pch = 19, col = colors)

# 添加从(0, 0)到(1, 1)的直线
abline(a=0, b=1, col="red") # a是截距，b是斜率

# 为每个点添加标签，使用最小的字体尺寸
text(data$NR, data$R, labels=data$Key, pos=4, cex=0.6)



pdf("G:\\TLSdata\\TLS_total_data\\ICB_varify\\GSE200996_scRNA\\Augur_plot.pdf",width = 4.5,height = 8) # Create a new pdf device
# 绘制散点图并为每个点使用分配的颜色
plot(data$R, data$delta_AUC, main = "Augur",
     xlab = "R", ylab = "delta_Augur_score", pch = 19, col = colors)


# 添加从(0, 0)到(1, 1)的直线
abline(a=0, b=0, col="red") # a是截距，b是斜率

# 为每个点添加标签，使用最小的字体尺寸
text(data$R, data$delta_AUC, labels=data$Key, pos=4, cex=0.6)
dev.off() # Close the pdf device





# 使用ggplot绘制散点图并添加文本
ggplot(data, aes(x = NR, y = R)) + 
  geom_point(aes(colour = colors)) +  # 绘制点
  geom_text_repel(aes(label = Key), box.padding = 0.35, point.padding = 0.5, max.overlaps = Inf) + 
  geom_abline(intercept = 0, slope = 1, col = "black") +  # 添加y=x的直线
  theme_minimal() +
  ggtitle("散点图示例")

# 添加调整后的图例
#legend("bottom", legend=names(rownames_colors), fill=rownames_colors, title="点的标识",
#       cex=0.8, ncol=3, x.intersp=0.5, y.intersp=0.5, inset=c(0,-0.1))
