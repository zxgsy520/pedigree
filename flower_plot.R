library(plotrix)

## read group.txt.type.xls file
args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
image_prefix = args[2]
cw = getwd()  #get current work directory
data = read.table(input_file,header=T,row.names=1,stringsAsFactors = F)
sample_id = colnames(data)
temp =colSums(data[c(3,5),])
samples_unique = c(temp)
#samples_unique = c(t(data[3,]))
samples_Single = data['Single',][1,1] 
sample_id_longest = max(nchar(sample_id))

if(sample_id_longest<=10){
	sample_id_plot_cex = 1
}else if(10<sample_id_longest & sample_id_longest<=15) {
	sample_id_plot_cex = 0.8
}else if(sample_id_longest>15){
	sample_id_plot_cex = 0.5
}

if(length(sample_id)<=20){
	ellipse_col = c('#6181BD4E','#F348004E','#64A10E4E','#9300264E','#464E044E','#049a0b4E','#4E0C664E','#D000004E','#FF6C004E','#FF00FF4E','#c7475b4E','#00F5FF4E','#BDA5004E','#A5CFED4E','#f0301c4E','#2B8BC34E','#FDA1004E','#54adf54E','#CDD7E24E','#9295C14E')
}else{
	ellipse_col = c(rep('#2B8BC34E', length(sample_id)))
}

flower_plot <- function(sample, unique_num, total, start, a, b, r, ellipse_col, circle_col) {
    par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
    plot(c(0,10),c(0,10),type='n')
    n   <- length(sample)
    deg <- 360 / n
    res <- lapply(1:n, function(t){
        draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180), 
            y = 5 + sin((start + deg * (t - 1)) * pi / 180), 
            col = ellipse_col[t],
            border = ellipse_col[t],
            a = a, b = b, angle = deg * (t - 1))
        text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
            y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
            unique_num[t])
        
        if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
            text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
                y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
                sample[t],
                srt = deg * (t - 1) - start,
                adj = 1,
                cex = sample_id_plot_cex
                )
        } else {
            text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
                y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
                sample[t],
                srt = deg * (t - 1) + start,
                adj = 0,
                cex = sample_id_plot_cex
                )
        }            
    })
    draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
    text(x = 5, y = 5, paste('Single:', total))
}

png(paste(image_prefix,".png",sep=""), width = 1500, height = 1500, res = 200, units = 'px')
flower_plot(sample = sample_id, unique_num = samples_unique , total = samples_Single,
    start = 90, a = 0.5, b = 2, r = 1, ellipse_col = ellipse_col, circle_col = 'white')
dev.off()
pdf(paste(image_prefix,".pdf",sep=""))
flower_plot(sample = sample_id, unique_num = samples_unique , total = samples_Single, 
    start = 90, a = 0.5, b = 2, r = 1, ellipse_col = ellipse_col, circle_col = 'white')
dev.off()
