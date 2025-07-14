library(hexSticker)


# make hex sticker ####
imgurl <- "~/Documents/earthstar.png"

sticker(subplot=imgurl, package = "specifiR",filename = "~/Desktop/hex_earthstar.png",
        s_width = .65,h_fill = "#948260",
        s_height = .65,
        s_x = 1,
        s_y = 1,
        p_y = 1.5,
        p_x = 1,
        p_size = 18,
        p_fontface = 'bold',
        p_color = "black",
        white_around_sticker = TRUE)



